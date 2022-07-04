#!/usr/bin/env python3
# coding: utf-8

# Split GBK
from Bio import SeqIO
import sys

for rec in SeqIO.parse(sys.argv[1], "genbank"):
    SeqIO.write([rec], open(rec.id + ".gbk", "w"), "genbank")
