#!/usr/bin/env python

'''This script takes the nextclade S gene protein translation file and counts the number of Xs and stop codons in the sequence'''

import sys

argvs = sys.argv

f = open(argvs[1])

max_N_prop = 0.01

seq_d = {}
for line in f:
  line = line.strip()
  if '>' in line:
    name = line.replace(">","")
    seq_d[name] = ""
  else:
    seq_d[name] += line

total_len = len(seq_d[name])

print("name\tlength\tcount_N\tcount_stop")
for name in seq_d:
  seq = seq_d[name]
  N_count = seq.count("X")
  total_len = len(seq)
  stop_count = seq.count("*")	
  print("\t".join([name,str(total_len),str(N_count),str(stop_count)]))
