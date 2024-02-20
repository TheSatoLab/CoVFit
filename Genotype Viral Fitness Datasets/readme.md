To generate the input data for CovFit training, running the make_input.sh bash script will perform the following operations:

Run the Nextclade CLI on the nucleotide sequence fasta from GISAID to genrate protein translations and find mutations
Run a python script to make a long-format table of the mutations present for each sequence from the nextclade output
Run a python script to filter out sequences based on the percentage of Xs and stop codons occuring in the sequence
Run an R script to calculate the Re of each sequence and output it along with the sequence's metadata
Run a python script that combines the Re and metadata information with the corresponding protein sequence to be used finally as the input for CovFit training
Required input data are the GISAID sequences fasta and metadata tsv files.
