# CovFit CLI
The CovFit CLI allows you to calculate fitness and DMS predictions on any standard fasta file of SARS-CoV-2 spike protein sequences from a Linux or WSL2 terminal (x86_64).

### Installation
Download the `CovFit_cli.tar.gz` file from `http...` and extract the contents.
Once extracted, you'll find the CovFit_cli executable inside a folder of the same name.

### Running
To view the full command line syntax and options for the CovFit CLI, `cd` to the `CovFit_cli` directory and input:<br>
`./CovFit_cli --help`

An example of using CovFit if you have a fasta file named `my_file.fasta` in your `~/Documents/Covid` directory would be:<br>
`./CovFit_cli --input ~/Documents/Covid/my_file.fasta --outdir ~/Documents/Covid/output/`<br>

The output file will be named `CovFit_Predictions_Fold_X.tsv`, where `X` is the number of the model instance used. Note that CovFit CLI will overwrite files in the output directory with the same fold number. 

### Folds
The CovFit CLI will run using one of the five model instances. The default is model instance 0. You can specify a particular model instance with the `-f, --fold' argument and an integer from 0 to 4.   

### DMS
By default, CovFit CLI will return predictions for viral fitness. Specifying the `-d, --dms` flag will additionally output the prediction scores for the 1,548 DMS tasks the model was trained on. 

### Batch
The batch size argument, set with `-b, --batch`, controls how many sequences are processed per batch. A high batch number may increase inference speed, but will increase memory consumption. The default batch size is 4.

### CPU vs GPU
By default, CovFit CLI will perform calculations on the CPU. Performing the calculations on a GPU using the `-g, --gpu` argument may significantly increase inference speed. An Nvidia GPU and CUDA installation are required. 

## Example with full options:
The following command:<br>
`./CovFit_cli --input my_file.fasta -outdir ~/Documents/Covid/output --fold 3 --dms --batch 16 --gpu` <br>
will perform inference on the GPU with 16 sequences per batch and output `~/Documents/Covid/output/CovFit_Predictions_Fold_3.tsv` where the first column is the mean fitness score accross countries, the next 1,548 columns are the DMS results, and the final 17 columns are the fitness results for each country.

To separate just the DMS results in a separate file, you can for instance use: `cut -f 1,3-1550 CovFit_Predictions_fold_3.tsv > just_DMS.tsv`, or for just the fitness fields: `cut --complement -f 3-1550 CovFit_Predictions_fold_3.tsv > just_fitness.tsv` 

