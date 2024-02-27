#!/usr/bin/env python3
# coding: utf-8

from datasets import Dataset
import pandas as pd
import numpy as np
import torch
from torch.optim import AdamW
from transformers import AutoTokenizer, DataCollatorForLanguageModeling, EsmForMaskedLM, get_cosine_schedule_with_warmup, Trainer, TrainingArguments 
from Bio import SeqIO as sio
import random
import math
import datetime


### Global Parameters

version = 'ESM2_650M_DA_combo_e30_lrsWarmupCosine'
model_type = 'facebook/esm2_t33_650M_UR50D'
max_length = 1000
batch_size = 5
warmup_epochs = 1
num_epochs = 30
randseed = 13


### Setting random seeds for reproducibility

torch.backends.cudnn.deterministic = True
random.seed(randseed)
torch.manual_seed(randseed)
torch.cuda.manual_seed(randseed)
np.random.seed(randseed)


### Data Prep

df_cv = pd.read_csv('df_coronaviradae_new.csv') # coronaviradae sequences
df_sc2 = pd.DataFrame.from_dict({i.id:str(i.seq) for i in sio.parse('legacy_seq_20220831_db99.fasta', 'fasta')},orient='index',columns=['seq']) # SARS-CoV-2 sequences
df_uniref = pd.read_table('../uniref/uniref_subset_seq.tsv') # uniref50 sequences


# Combine cornoviradae and SARS-CoV-2 sequences to make training set
X_train = df_cv.seq.tolist()+df_sc2.seq.tolist()
train_df = pd.DataFrame(X_train)
train_df.columns=['seq']


# Load main dataset and split data by date to make an evaluation dataset of "future" sequences
full_df = pd.read_table('metadata.representative.all_countries.with_date.with_aligned_seq.with_immune_escape_score.txt')
input_df = full_df[['relative_Re', 'country', 'seq', 'date.first','relative_Re_log', 'clade']]
cutoff = pd.to_datetime('2022-09-01')
input_df['date.first'] = pd.to_datetime(input_df['date.first'])
recent_df = input_df[input_df['date.first']>=cutoff]


# Input data into dataset.Dataset to prevent incompatibility with slow tokenizer and data_collator 
train_dataset = Dataset.from_pandas(train_df[['seq']]) 
recent_dataset = Dataset.from_pandas(recent_df[['seq']])
test_dataset = Dataset.from_pandas(df_uniref[['seq']])


# Tokenize the datasets
tokenizer = AutoTokenizer.from_pretrained(model_type)

def tokenizer_for_map(input): #Tokenizer and params including special_tokens_mask required for MLM
    return tokenizer(
        input['seq'],
        padding='max_length',
        truncation=True,
        max_length=1000,
        return_special_tokens_mask=True
    )
  
column_names = train_dataset.column_names #This will be the names of all the old columns, to then be deleted after the new tokenized columns are added.

train_dataset = train_dataset.map( 
    tokenizer_for_map,
    batched=True,
    num_proc=8,
    remove_columns=column_names,
)

recent_dataset = recent_dataset.map(
    tokenizer_for_map,
    batched=True,
    num_proc=8,
    remove_columns=column_names,
)

test_dataset = test_dataset.map(
    tokenizer_for_map,
    batched=True,
    num_proc=8,
    remove_columns=column_names,
)


### Training

data_collator = DataCollatorForLanguageModeling(tokenizer=tokenizer,return_tensors='pt',mlm_probability=0.15) # Provides random masking and returns tensors during training per-batch
model = EsmForMaskedLM.from_pretrained(model_type) # Model with MaskedLM layer on top to provide predictions of masked sequences

params = filter(lambda x: x.requires_grad, model.parameters())
optimizer = torch.optim.AdamW(params, lr=2e-5)
warmup_steps = math.ceil((train_dataset.num_rows/batch_size)*warmup_epochs)
training_steps = math.ceil((train_dataset.num_rows/batch_size)*num_epochs)
scheduler = get_cosine_schedule_with_warmup(optimizer, num_warmup_steps=warmup_steps, num_training_steps=training_steps)


training_args = TrainingArguments(
    output_dir=f"{version}/trainer",
    overwrite_output_dir=True,
    num_train_epochs=num_epochs,
    per_device_train_batch_size=batch_size,
    evaluation_strategy='no',
    do_eval=False,
    save_total_limit=1,
    seed=randseed,
    dataloader_num_workers=4,
)

trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=train_dataset,
    data_collator=data_collator,
    optimizers=[optimizer,scheduler]
)

trainer.train()

trainer.save_model(f"{version}/model")


### Evaluation on uniref50

model = EsmForMaskedLM.from_pretrained(model_type)

evaluator = Trainer(
  model=model,
  data_collator=data_collator,
  eval_dataset=test_dataset,
  )

eval_results_old = evaluator.evaluate()


model = EsmForMaskedLM.from_pretrained(f"{version}/model")

evaluator = Trainer(
  model=model,
  data_collator=data_collator,
  eval_dataset=test_dataset,
  )

eval_results_new = evaluator.evaluate()


import math
print('ESM2 Uniref Data:')
print(f"Original Model's Perplexity: {math.exp(eval_results_old['eval_loss']):.4f}")
print(f"Adapted Model's Perplexity: {math.exp(eval_results_new['eval_loss']):.4f}")


### Evaluation on "future" SARS-CoV-2 sequences

model = EsmForMaskedLM.from_pretrained(model_type)

evaluator = Trainer(
  model=model,
  data_collator=data_collator,
  eval_dataset=recent_dataset,
  )

recent_results_old = evaluator.evaluate()


model = EsmForMaskedLM.from_pretrained(f"{version}/model")

evaluator = Trainer(
  model=model,
  data_collator=data_collator,
  eval_dataset=recent_dataset,
  )

recent_results_new = evaluator.evaluate()


print('Future SARS-CoV-2 Data:')
print(f"Original Model's Perplexity: {math.exp(recent_results_old['eval_loss']):.4f}")
print(f"Adapted Model's Perplexity: {math.exp(recent_results_new['eval_loss']):.4f}")
