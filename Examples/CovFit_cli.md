# CovFit CLI
The CovFit CLI allows you to calculate fitness and DMS predictions on any standard fasta file of SARS-CoV-2 spike protein sequences from a Linux or WSL2 terminal.

### Installation
Download the `CovFit_cli.tar.gz` file from `http...` and extract the contents.
Once extracted, you'll find the CovFit_cli executable inside a folder of the same name.

### Running
To view the command line syntax and options for the CovFit CLI, `cd` to the `CovFit_cli` directory and input:
`./CovFit_cli --help`

An example of using CovFit if you have a fasta file named `my_file.fasta` in your `~/Documents/Covid` directory would be:
`./CovFit_cli --input ~/Documents/Covid/my_file.fasta --outdir ~/Documents/Covid/output/`


