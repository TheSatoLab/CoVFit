#!/usr/bin/env python

'''This script takes the nextclade S gene protein translation file and counts the number of Xs and stop codons in the sequence'''

import sys
from adam_utils.bio import fasta

argvs = sys.argv

f = open(argvs[1])

#max_N_prop = 0.01 # This setting and filtering is actually in the R script

seq_d = fasta.to_dict(f,unique=True)

#total_len = len(seq_d[name])

print("name\tlength\tcount_N\tcount_stop")
for name in seq_d:
  seq = seq_d[name]
  N_count = seq.count("X")
  total_len = len(seq)
  stop_count = seq.count("*")	
  print("\t".join([name,str(total_len),str(N_count),str(stop_count)]))
