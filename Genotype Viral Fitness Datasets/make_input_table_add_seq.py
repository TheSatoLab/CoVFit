#!/usr/bin/env python

import sys
import pandas as pd
from Bio import SeqIO as sio

def to_dict(path=None, file_format='fasta', robust=False, unique=False):
    if path is None:
        print('Path to file required')
        return
    if robust:
        print('Forced unique sequence names')
        # `record.description` を使用してヘッダー全体を取得
        res = {f'{i}_'+record.description: str(record.seq) for i, record in enumerate(sio.parse(path, file_format))}
    elif unique:
        print('Removed duplicate name occurrences without checking sequences')
        names = set()
        res = {}
        for record in sio.parse(path, file_format):
            # ここでも `record.description` を使用
            if record.description not in names:
                names.add(record.description)
                res[record.description] = str(record.seq)
    else:
        # ヘッダー全体をキーとして使用
        res = {record.description: str(record.seq) for record in sio.parse(path, file_format)}

    return res
argvs = sys.argv

metadata_f = open(argvs[1])
fasta_f = open(argvs[2])

header = metadata_f.__next__().strip()

name_l = []
metadata_l = []
for line in metadata_f:
  line = line.strip().split("\t")
  name = line[3]
  name_l.append(name)
  metadata_l.append(line)


name_s = set(name_l)

seq_d = to_dict(path=fasta_f, unique=True)

print(header + "\tseq")
for i in range(len(name_l)):
  name = name_l[i]
  line = metadata_l[i]
  seq = seq_d[name] #.replace("-","")
  line.append(seq)
  print("\t".join(line))
