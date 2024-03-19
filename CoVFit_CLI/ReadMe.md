# CovFit CLI
The CovFit CLI allows you to calculate fitness and DMS predictions on any standard fasta file of SARS-CoV-2 spike protein sequences from a Linux or Windows WSL2 terminal (x86_64). The CLI was packaged on Ubuntu 20.04 and has been tested on 18.04 and 22.04.

## Installation
Download the `covfit_cli.tar.gz` file from Zenodo at `http...` and extract the contents. You may be able to extract the tar.gz file by double-clicking on it in your file explorer. Alternatively, it can be extracted from the command line with: <br>
`tar -xzvf path/to/folder/covfit_cli.tar.gz` <br>
Once extracted, you'll find the covfit_cli executable inside a folder of the same name. No specific installation procedure is required beyond extracting the files.

## Running
To view the full command line syntax and options for the CoVFit CLI, `cd` to the `CoVFit_cli` directory and input:<br>
`./covfit_cli --help`
#### Test
A sample fasta file is included in the CoVFit directory for testing purposes. Run CoVFit on the sample sequences with the following simple command: `./covfit_cli -i covfit_samples.fasta`, which will output `CoVFit_Predictions_Fold_0.tsv`. The CoVFit directory also includes `covfit_samples_nextclade.tsv` which includes nextclade metadata for the sample sequences.  
#### Example
If you have a fasta file named `my_file.fasta` in your `~/Documents/covid19/` directory would be:<br>
`./covfit_cli --input ~/Documents/covid19/my_file.fasta --outdir ~/Documents/covid19/output/`<br>
 The output file will be named `CoVFit_Predictions_Fold_X.tsv`, where `X` is the number of the model instance used. The file will be placed in the location designated by `-o, --outdir`, in this example `~/Documents/covid19/output/`. If the output directory does not exist, it will be created. <br>

A fasta file with 100 sequences can generally be processed within a couple minutes, depending on system performance. **Note that CoVFit CLI will overwrite files in the output directory with the same fold number**. 

### Folds
The CoVFit CLI will run using one of the five model instances. The default is model instance 0. You can specify a particular model instance with the `-f, --fold' argument and an integer from 0 to 4.   

### DMS
By default, CoVFit CLI will return predictions for viral fitness. Specifying the `-d, --dms` flag will additionally output the prediction scores for the 1,548 DMS tasks the model was trained on. 

### Batch
The batch size argument, set with `-b, --batch`, controls how many sequences are processed per batch. A high batch number may increase inference speed, but will increase memory consumption. The default batch size is 4.

### CPU vs GPU
By default, CoVFit CLI will perform calculations on the CPU. Performing the calculations on a GPU using the `-g, --gpu` argument may significantly increase inference speed. An Nvidia GPU and CUDA installation are required. 

## Example with full options:
The following command:<br>
`./covfit_cli --input my_file.fasta -outdir ~/Documents/covid19/output/ --fold 3 --dms --batch 16 --gpu` <br>
will perform inference on the GPU with 16 sequences per batch and output the file `~/Documents/covid19/output/CoVFit_Predictions_Fold_3.tsv` where the first column is the mean fitness score accross countries, the next 1,548 columns are the DMS results, and the final 17 columns are the fitness results for each country.

To separate just the DMS results in a separate file, you can for instance use: `cut -f 1,3-1550 CoVFit_Predictions_fold_3.tsv > just_DMS.tsv`, or for just the fitness fields: `cut --complement -f 3-1550 CoVFit_Predictions_fold_3.tsv > just_fitness.tsv` 

