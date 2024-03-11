# CovFit CLI
The CovFit CLI allows you to calculate fitness and DMS predictions on any standard fasta file of SARS-CoV-2 spike protein sequences from a Linux or WSL2 terminal.

### Installation
Download the `CovFit_cli.tar.gz` file from `http...` and extract the contents.
Once extracted, you'll find the CovFit_cli executable inside a folder of the same name.

### Running
To view the full command line syntax and options for the CovFit CLI, `cd` to the `CovFit_cli` directory and input:<br>
`./CovFit_cli --help`

An example of using CovFit if you have a fasta file named `my_file.fasta` in your `~/Documents/Covid` directory would be:<br>
`./CovFit_cli --input ~/Documents/Covid/my_file.fasta --outdir ~/Documents/Covid/output/`<br>
The output file will be named CovFit_Predictions_Fold_X.tsv, where X is the number of the model instance used. 

### Folds
The CovFit CLI will run using one of the five model instances. The default is model instance 0. You can specify a particular model instance with the `-f, --fold' argument and an integer from 0 to 4.   

### DMS
By default, CovFit CLI will return predictions for viral fitness. Specifying the `-d, --dms` flag will additionally output the prediction scores for the 1,548 DMS tasks the model was trained on. 

### Batch
The batch size argument, set with `-b, --batch`, controls how many sequences are processed per batch. A high batch number may increase inference speed, but will increase memory consumption. The default batch size is 4.

### CPU vs GPU
By default, CovFit CLI will perform calculations on the CPU. Performing the calculations on a GPU using the `-g, --gpu` argument may significantly increase inference speed. An Nvidia GPU and CUDA installation are required. 

## Example with full options:
`./CovFit_cli --input my_file.fasta -outdir ~/Documents/Covid/output --fold 3 --dms --batch 16 --gpu` 

