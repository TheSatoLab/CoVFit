#!/usr/bin/env python
'''Takes the nextclade output file and filters it down to name and mutations'''
import sys

argvs = sys.argv

f = open(argvs[1])

f.__next__()

print("seqName\tmut")
for line in f:
  line = line.strip().split("\t")
  name = line[1]
  s,d,i = line[29:32]
  for m in [s,d,i]:
    m_l = m.split(",")
    for m_i in m_l:
      if m_i != "":
        print("\t".join([name,m_i]))
