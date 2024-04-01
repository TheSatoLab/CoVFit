#!/usr/bin/env python

import sys
import random
import glob
import re
import numpy as np
from tqdm import tqdm
import pandas as pd
import torch
from torch.utils.data import DataLoader, Dataset, TensorDataset
from transformers import AutoTokenizer, AutoConfig, EsmModel, EsmForSequenceClassification, AutoModel, AdamW
import os
import torch.nn as nn
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
from Bio import SeqIO
from sklearn.preprocessing import MinMaxScaler
import bitsandbytes as bnb

argvs = sys.argv

#hyper parameters###############################
max_length = 1024 # max token length
cutoff = pd.to_datetime("2022-09-01") # cutoff date
min_country = 300 # Countries with <300 genotype-fitness data were removed.
date_thresh = pd.to_datetime("2022-01-01") # for loss weighting by variant emergence date
target_strain = "BA.2"
negative_clade_list = ["22E","22F","23A","23B","23C","23D","23E","23F","BA.2.86","recombinant"]

#inputs###############################
dir="/work/ge17/e17000/variant_fitness/esm_multitask/"
MODEL_NAME = "facebook/esm2_t33_650M_UR50D"
DA_model_name = dir + "DA_model/model_uniESM2_650M_DA_combo_e30_lrsWarmupCosine"

fold_id = int(argvs[1])
output_prefix = dir + "output/" + argvs[2]

full_df_name = dir + "input_230829/metadata.representative.all_countries.with_date.v2.with_seq.txt"
antibody_escape_f_name = dir + "input_230829/escape_data_mutation.csv"
dms_fasta_f_name = dir + "input_230829/nextclade.peptide.S_rename.fasta"

### Setting random seeds for reproducibility
randseed = 13
torch.backends.cudnn.deterministic = True
random.seed(randseed)
torch.manual_seed(randseed)
torch.cuda.manual_seed(randseed)
np.random.seed(randseed)



#functions###############################
#reading fasta
def fasta_to_dataframe(fasta_file_name):
    descriptions = []
    sequences = []
    for record in SeqIO.parse(fasta_file_name, "fasta"):
        descriptions.append(record.description)
        sequences.append(str(record.seq))

    seq_df = pd.DataFrame({
        'target': descriptions,
        'seq': sequences
    })
    return seq_df


#inserting mutations to generate mutant S proteins for DMS data
def mutate_sequence(row):
    site = int(row["position"])
    seq = row["seq"]
    mut = row["mutant"]
    return seq[:site-1] + mut + seq[site:]


#scaling DMS data
def scale_values(x):
    # 95パーセンタイルを取得
    max_value = np.percentile(x, 95)
    min_value = x.min()
    return (x - min_value) / (max_value - min_value)


#generating task_info dictionary
def encode_categorical_variables(categories):
    unique_categories = sorted(set(categories))
    category_to_number = {category: i for i, category in enumerate(unique_categories)}
    numerical_categories = [category_to_number[category] for category in categories]
    return numerical_categories, category_to_number


#splitting training data
from sklearn.model_selection import StratifiedKFold

def stratified_kfold_data_split(df, k=5, fold_id=0, seed=randseed):
    y = df['split_class']
    kf = StratifiedKFold(n_splits=k, shuffle=True, random_state=seed)
    splits = list(kf.split(df, y))

    if fold_id < 0 or fold_id >= k:
        raise ValueError("fold_id should be between 0 and n_splits-1")

    train_val_index, test_index = splits[fold_id]
    train_val = df.iloc[train_val_index].reset_index(drop=True)
    test = df.iloc[test_index].reset_index(drop=True)

    y_train_val = train_val['split_class']
    kf_train_val = StratifiedKFold(n_splits=4, shuffle=True, random_state=seed)
    splits_train_val = list(kf_train_val.split(train_val, y_train_val))
    train_index, val_index = splits_train_val[0]
    train = train_val.iloc[train_index].reset_index(drop=True)
    val = train_val.iloc[val_index].reset_index(drop=True)

    return train, val, test




##making data chunks, composed of a tokenized input seq and outcome variables for multiple tasks (up to 10 tasks)
from joblib import Parallel, delayed

#tokenize
def parallel_tokenizer_chunk(seq_chunk, tokenizer, max_length):
    return tokenizer(
        seq_chunk,
        max_length=max_length,
        padding="max_length",
        truncation=True,
        return_attention_mask=True,
        return_tensors="pt"
    )

#non-overlapping sampling
import random
def sample_and_replace(outcome_l, weight_l, idx_l, n):
    index_l = [i for i, value in enumerate(outcome_l) if not np.isnan(value)]
    result_outcome_l = []
    result_weight_l = []
    result_idx_l = []
    while index_l:
        sampled_index_l = random.sample(index_l, min(n, len(index_l)))

        sampled_outcome_l = outcome_l.copy()
        for i in range(len(sampled_outcome_l)):
            if i not in sampled_index_l:
                sampled_outcome_l[i] = np.nan
        result_outcome_l.append(sampled_outcome_l)

        sampled_weight_l = weight_l.copy()
        for i in range(len(sampled_weight_l)):
            if i not in sampled_index_l:
                sampled_weight_l[i] = np.nan
        result_weight_l.append(sampled_weight_l)

        sampled_idx_l = idx_l.copy()
        for i in range(len(sampled_idx_l)):
            if i not in sampled_index_l:
                sampled_idx_l[i] = np.nan
        result_idx_l.append(sampled_idx_l)
        index_l = [i for i in index_l if i not in sampled_index_l]

    return result_outcome_l, result_weight_l, result_idx_l

#making data chunks
def preprocess_data(df, n_targets, max_length, tokenizer, n_samples = 10):
    nans = [np.nan] * n_targets

    unique_seq_ids = df['seq_id'].unique()
    outcome_sum_d = {seq_id: nans.copy() for seq_id in unique_seq_ids}
    weight_sum_d = {seq_id: nans.copy() for seq_id in unique_seq_ids}
    idx_sum_d = {seq_id: nans.copy() for seq_id in unique_seq_ids}

    for index, row in df.iterrows():
        seq_id, task_id, outcome, weight = row['seq_id'], row['task_id'], row['outcome'], row['weight']
        outcome_sum_d[seq_id][task_id] = outcome
        weight_sum_d[seq_id][task_id] = weight
        idx_sum_d[seq_id][task_id] = index

    unique_seq_df = df.drop_duplicates(subset='seq_id')[['seq_id','seq']]

    # Split sequences into chunks for parallel processing
    n_chunks = 10
    seq_chunks = [list(chunk) for chunk in np.array_split(unique_seq_df["seq"].tolist(), n_chunks)]

    encodings = Parallel(n_jobs=-1)(delayed(parallel_tokenizer_chunk)(seq_chunk, tokenizer, max_length) for seq_chunk in seq_chunks)

    # Combine the results from all chunks
    input_ids = torch.cat([encoding["input_ids"] for encoding in encodings], dim=0)
    attention_masks = torch.cat([encoding["attention_mask"] for encoding in encodings], dim=0)

    outcomes = [outcome_sum_d[seq_id] for seq_id in unique_seq_ids]
    weights = [weight_sum_d[seq_id] for seq_id in unique_seq_ids]
    idxs = [idx_sum_d[seq_id] for seq_id in unique_seq_ids]

    input_ids_final = []
    attention_masks_final = []
    outcomes_final = []
    weights_final = []
    idxs_final = []

    for i in range(len(input_ids)):
        input_id = input_ids[i]
        attention_mask = attention_masks[i]
        outcome = outcomes[i]
        weight = weights[i]
        idx = idxs[i]
        outcome2, weight2, idx2 = sample_and_replace(outcome, weight, idx, n_samples)
        input_id2 = [input_id] * len(outcome2)
        attention_mask2 = [attention_mask] * len(outcome2)
        input_ids_final.extend(input_id2)
        attention_masks_final.extend(attention_mask2)
        outcomes_final.extend(outcome2)
        weights_final.extend(weight2)
        idxs_final.extend(idx2)

    outcomes_final = torch.tensor(outcomes_final, dtype=torch.float32)
    weights_final = torch.tensor(weights_final, dtype=torch.float32)
    idxs_final = torch.tensor(idxs_final, dtype=torch.float32)
    return input_ids_final, attention_masks_final, outcomes_final, weights_final, idxs_final



#checking trainable parameters (for LoRA)
def print_trainable_parameters(model):
    """
    Prints the number of trainable parameters in the model.
    """
    trainable_params = 0
    all_param = 0
    for _, param in model.named_parameters():
        all_param += param.numel()
        if param.requires_grad:
            trainable_params += param.numel()
    print(
        f"trainable params: {trainable_params} || all params: {all_param} || trainable%: {100 * trainable_params / all_param:.2f}"
    )




#reading dms backbone sequence###############################

dms_seq_df = fasta_to_dataframe(dms_fasta_f_name)

#reading, processing, filtering antibody escape DMS###############################
antibody_escape_df = pd.read_csv(antibody_escape_f_name)

negative_source_l = ['SARS convalescents','WT-engineered'] # filtering by source
antibody_escape_df = antibody_escape_df[~antibody_escape_df['source'].isin(negative_source_l)]
antibody_escape_df = antibody_escape_df[antibody_escape_df["IC50"]<10] #filtering by IC50

antibody_escape_df = antibody_escape_df.rename(columns={'target_virus': 'target', 'mutation':'mutant', 'site':'position'})

antibody_escape_df = antibody_escape_df[antibody_escape_df["target"].isin([target_strain])] #filtering by target variant

antibody_escape_df.loc[:,'group'] = antibody_escape_df["target"] + "_" + antibody_escape_df["group"] #group = epitope group; group = task_group (used for balancing los weights)
antibody_escape_df.loc[:,'data_group'] = antibody_escape_df["target"] + "_" + antibody_escape_df["condition"] #condition = mAb type; data_group = task_name

antibody_escape_df.loc[:,'mut_escape_w'] = antibody_escape_df["mut_escape"] * antibody_escape_df["neg_log_IC50"] #weighting by IC50

antibody_escape_df.loc[:,"outcome_pre"] = np.log10(antibody_escape_df["mut_escape_w"]+1) # log transformation with pusedo count 1


antibody_escape_df.loc[:,'outcome'] = antibody_escape_df.groupby('group')['outcome_pre'].transform(scale_values) # scaling so that 0 and 95 percentile fell within the range 0–1
antibody_escape_df.loc[:,'outcome'] = antibody_escape_df['outcome'].clip(0, 1) #clipping >1 value to 1

#making S mutant sequences
antibody_escape_df = antibody_escape_df.merge(dms_seq_df, on='target')
antibody_escape_df.loc[:,"seq"] = antibody_escape_df.apply(mutate_sequence, axis=1)

antibody_escape_df_selected = antibody_escape_df[['group','data_group','outcome','seq']]


#reading, filtering, processing genotype-fitness data###############################
full_df = pd.read_table(full_df_name)
full_df = full_df.groupby('country').filter(lambda x: len(x) >= min_country) #removing countires with fewer data

full_df.loc[:,'clade'] = np.where(full_df["Nextclade_pango"].isin(["BA.2.86","BA.2.86.1","JN.1","JN.2"]), "BA.2.86", full_df['clade']) #renaming BA.2.86

full_df.loc[:,'relative_Re_log'] = np.log(full_df['relative_Re']) #log transformation

#scaling so that the 0.1 percentile and 99.9 percentile points fall between 0 and 1
q1 = full_df['relative_Re_log'].quantile(0.001)
q99 = full_df['relative_Re_log'].quantile(0.999)

full_df.loc[:, 'relative_Re_norm'] = (full_df['relative_Re_log'] - q1) / (q99 - q1)



input_df = full_df[['relative_Re', 'country', 'seq', 'date.first','relative_Re_norm', 'clade']]

input_df.loc[:,"group"] = "fitness" # group = task_group (used for balancing los weights)
input_df.loc[:,"data_group"] = "fitness_" + input_df["country"] # data_group = task_name

input_df_selected = input_df[["group","data_group","relative_Re_norm","seq",'date.first','clade','country']]
input_df_selected.columns = ["group","data_group","outcome","seq",'date.first','clade','country']

input_df_selected.loc[:,'date.first']=pd.to_datetime(input_df_selected['date.first'])


#combining data###############################
cat_df = pd.concat([input_df_selected, antibody_escape_df_selected])
cat_df = cat_df.sample(frac=1,random_state=randseed)
cat_df = cat_df.set_axis(range(len(cat_df)))

#removing genotypes with >5 Xs or >30 hyphens###############################
count_x = [c.count("X") for c in cat_df["seq"]]
count_hyphen = [c.count("-") for c in cat_df["seq"]]

cat_df.loc[:,"count_x"] = count_x
cat_df.loc[:,"count_hyphen"] = count_hyphen

cat_df = cat_df[cat_df["count_x"]<=5]
cat_df = cat_df[cat_df["count_hyphen"]<=30]


#making task_id_infos and seq_id_info###############################
task_ids,task_id_infos = encode_categorical_variables(cat_df["data_group"].tolist())

cat_df.loc[:,"task_id"] = task_ids


seq_ids,seq_id_infos = encode_categorical_variables(cat_df["seq"].tolist())

cat_df.loc[:,"seq_id"] = seq_ids


#balancing loss weights between fitness and DMS tasks###############################
cat_df.loc[:,'major_group'] = np.where(cat_df['group'] == 'fitness', 'fitness','immune')

freq = cat_df['major_group'].value_counts()
freq_dict = freq.to_dict()

num_fitness_past = len(cat_df[cat_df['date.first']<cutoff])


freq_dict['fitness'] = num_fitness_past

weight_dict = {}
for k,v in freq_dict.items():
  w = 1 / v
  weight_dict[k] = w

w_max = max(weight_dict.values())

for k,v in weight_dict.items():
  v = v / w_max
  weight_dict[k] = v

cat_df.loc[:,'major_group2'] = np.where((cat_df['major_group'] == 'fitness') & (cat_df['date.first'] >= date_thresh), 'fitness_recent',
                                     np.where((cat_df['major_group'] == 'fitness') & (cat_df['date.first'] < date_thresh), 'fitness',cat_df['major_group']))

weight_dict["fitness"] = weight_dict["fitness"] * 1 # In the past-future experiments, weights only for variants emrged later than 2022-01-01 were doubled.
weight_dict["fitness_recent"] = weight_dict["fitness"] * 2 # doubling weights for recent variants

cat_df.loc[:,"weight"] = [weight_dict[c] for c in cat_df['major_group2']]


#data splitting###############################
#split class = (mAb type for DMS data) or (the combination of genotypes and countries for genotype-fitness data)
cat_df.loc[:,"split_class"] =  np.where((cat_df['major_group'] == 'fitness'), cat_df['data_group'] + "_" + cat_df['clade'], cat_df['data_group'])

#preserving future
future_df = cat_df[cat_df['date.first']>=cutoff]
cat_df = cat_df[(cat_df['date.first']<cutoff) | (cat_df['date.first'].isna())]
cat_df = cat_df[~cat_df["clade"].isin(negative_clade_list)] #removing future variants of interest from the past dataset

#data split
cat_df = cat_df.sample(frac=1,random_state=randseed)
cat_df = cat_df.set_axis(range(len(cat_df)))

cat_df_train, cat_df_val, cat_df_test = stratified_kfold_data_split(cat_df, fold_id=fold_id)
future_df = future_df.reset_index(drop=True)

#preparing dataset###############################
from transformers import AutoTokenizer, Trainer, TrainingArguments, DefaultDataCollator

tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)

n_targets = len(task_id_infos)

input_ids_train, attention_masks_train, outcomes_train, weights_train, idxs_train = preprocess_data(cat_df_train, n_targets, max_length, tokenizer)
input_ids_val, attention_masks_val, outcomes_val, weights_val, idxs_val = preprocess_data(cat_df_val, n_targets, max_length, tokenizer)
input_ids_test, attention_masks_test, outcomes_test, weights_test, idxs_test = preprocess_data(cat_df_test, n_targets, max_length, tokenizer)
input_ids_future, attention_masks_future, outcomes_future, weights_future, idxs_future = preprocess_data(future_df, n_targets, max_length, tokenizer, n_samples = 17)

class TextDataset(Dataset):
    def __init__(self, input_ids, attention_masks, labels, weights):
        self.input_ids = input_ids
        self.attention_masks = attention_masks
        self.labels = labels
        self.weights = weights

    def __getitem__(self, idx):
        item = {
            'input_ids': self.input_ids[idx],
            'attention_masks': self.attention_masks[idx],
            'labels': self.labels[idx],
            'weights': self.weights[idx]
        }
        return item

    def __len__(self):
        return len(self.labels)


dataset_train = TextDataset(input_ids_train, attention_masks_train, outcomes_train, weights_train)
dataset_val = TextDataset(input_ids_val, attention_masks_val, outcomes_val, weights_val)
dataset_test = TextDataset(input_ids_test, attention_masks_test, outcomes_test, weights_test)
dataset_future = TextDataset(input_ids_future, attention_masks_future, outcomes_future, weights_future)

#model construction###############################
from transformers import EsmModel, EsmConfig, PreTrainedModel
from transformers.modeling_outputs import SequenceClassifierOutput
import torch
import torch.nn as nn

class EsmForRegression(PreTrainedModel):
    config_class = EsmConfig
    #setting task-specific regression heads
    def __init__(self, config, n_targets, intermediate_dim=256, dropout_rate=0.1):
        super(EsmForRegression, self).__init__(config)
        self.esm = EsmModel(config)
        self.regressor = nn.Sequential(
            nn.Linear(config.hidden_size, intermediate_dim),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(intermediate_dim, n_targets)
        )
    #weighted_mse_loss
    def weighted_mse_loss(self, inputs, targets, weights):
        diff = inputs - targets
        diff_squared = diff ** 2
        weighted_diff_squared = diff_squared * weights
        loss = weighted_diff_squared.mean()
        return loss

    def forward(self, input_ids, attention_mask=None, labels=None, weights=None, **kwargs):
        outputs = self.esm(input_ids, attention_mask=attention_mask)
        pooled_output = outputs[0][:, 0]
        regression_output = self.regressor(pooled_output)

        if weights is None:
            weights = torch.ones_like(regression_output)

        #calculating loss only for tasks included in a data chunk
        loss = None
        if labels is not None:
            mask = ~torch.isnan(labels)
            masked_labels = labels[mask]
            masked_regression_output = regression_output[mask]
            masked_weights = weights[mask]
            loss = self.weighted_mse_loss(masked_regression_output, masked_labels, masked_weights)

        return SequenceClassifierOutput(
            loss=loss,
            logits=regression_output,
            hidden_states=outputs.hidden_states,
            attentions=outputs.attentions,
    )


esm_config = EsmConfig.from_pretrained(DA_model_name)
model = EsmForRegression(esm_config,n_targets)
model.esm = AutoModel.from_pretrained(DA_model_name)

#setting LoRA###############################
from peft import LoraConfig, get_peft_model

config = LoraConfig(
    task_type="SEQ_CLS",
    r=8,
    lora_alpha=16,
    target_modules=["key", "query", "value", "dense"],
    lora_dropout=0.05,
    bias="lora_only",
    modules_to_save=["regressor"]
)
lora_model = get_peft_model(model, config)
print_trainable_parameters(lora_model)


#training###############################
from scipy.stats import pearsonr, spearmanr
from transformers import Trainer, TrainingArguments

def compute_metrics(eval_pred):
    predictions = torch.tensor(eval_pred.predictions)
    labels = torch.tensor(eval_pred.label_ids)

    mask = ~torch.isnan(labels)
    masked_labels = labels[mask]
    masked_predictions = predictions[mask]
    pearson_corr = pearsonr(masked_predictions, masked_labels)[0]
    spearman_corr = spearmanr(masked_predictions, masked_labels)[0]

    return {
        'pearson': pearson_corr,
        'spearman': spearman_corr
    }


training_args = TrainingArguments(
    output_dir= output_prefix + "results",
    evaluation_strategy="epoch",
    learning_rate=2e-4,
    per_device_train_batch_size=4,
    per_device_eval_batch_size=16,
    gradient_accumulation_steps=2,
    fp16=True,
    load_best_model_at_end=True,
    remove_unused_columns=False,
    logging_steps=1,
    num_train_epochs=30,
    weight_decay=0.02,
    logging_dir="./logs",
    save_strategy="epoch",
    save_total_limit=1,
    push_to_hub=False,
)

trainer = Trainer(
    model=lora_model,
    args=training_args,
    train_dataset=dataset_train,
    eval_dataset=dataset_val,
    compute_metrics=compute_metrics,
)

trainer.train()

#save model###############################
model_f_name = output_prefix + "model.ckpt"
torch.save(lora_model.state_dict(), model_f_name)


#prediction###############################
predictions = trainer.predict(dataset_test)

mask = ~torch.isnan(idxs_test)

masked_outcomes_test = outcomes_test[mask]
masked_predictions = torch.tensor(predictions.predictions)[mask]
masked_idxs_test = idxs_test[mask].int()


predicted_df = pd.DataFrame({
    'predicted_outcome': masked_predictions.numpy()
})
predicted_df.index = masked_idxs_test.tolist()

predicted_df = pd.merge(predicted_df, cat_df_test, left_index=True, right_index=True)


#prediction for future variants###############################
predictions_future = trainer.predict(dataset_future)

mask_future = ~torch.isnan(idxs_future)

masked_outcomes_future = outcomes_future[mask_future]
masked_predictions_future = torch.tensor(predictions_future.predictions)[mask_future]
masked_idxs_future = idxs_future[mask_future].int()


predicted_df_future = pd.DataFrame({
    'predicted_outcome': masked_predictions_future.numpy()
})
predicted_df_future.index = masked_idxs_future.tolist()
predicted_df_future = pd.merge(predicted_df_future, future_df, left_index=True, right_index=True)
predicted_df_future_fitness = predicted_df_future[predicted_df_future["group"] == "fitness"]

#output###############################
out_test_f_name = output_prefix + "predicted_test.txt"
predicted_df.to_csv(out_test_f_name,header=True, index=False, sep='\t')

out_future_f_name = output_prefix + "predicted_future.txt"
predicted_df_future_fitness.to_csv(out_future_f_name,header=True, index=False, sep='\t')

