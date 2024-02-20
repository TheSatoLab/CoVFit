nextclade run -d sars-cov-2 --genes S --output-tsv=full.nextclade.tsv --output-translations='output_dir/gene_{gene}.translation.fasta' sequences.fasta

python make_mut_table.py full.nextclade.tsv > full.nextclade.mut_long.tsv

python count_X.py gene_S.translation.fasta > gene_S.translation.count_table.txt

Rscript make_input_all_region.3.R

python make_input_table_add_seq.py metadata.representative.all_countries.with_date.txt gene_S.translation.fasta >  metadata.representative.all_countries.with_date.with_aligned_seq.txt
