#!/usr/bin/env python

# impor required libraries
from Bio import SeqIO
import sys

# get filename from command line
filename = sys.argv[-1]
gbk_filename = filename

# open file connection
input_handle  = open(gbk_filename, "r")

# def function for stderr print
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

for seq_record in SeqIO.parse(input_handle, "genbank") :
    
    # all brank to check for errors later
    organism = None
    product  = None
    gene     = None
    
    for seq_feature in seq_record.features:

        # get organism
        if seq_feature.type.lower()=="source" :
            organism = seq_feature.qualifiers['organism'][0]

        # get product
        if seq_feature.type.lower()=="protein" :
            product = seq_feature.qualifiers['product'][0].replace(" ", "_")
        
        # get gene name if entry has only CDS definition
        if seq_feature.type.lower()=="cds" :
            gene = seq_feature.qualifiers['gene'][0]
        
        # get gene name if entry has only Gene definition
        if seq_feature.type.lower()=="gene" :
            gene = seq_feature.qualifiers['gene'][0]

        # save sequence info
        seq  = seq_record.seq
        acc  = seq_record.name
    
    # print
    if gene==None or product==None:
        eprint(f"An error has been found with entry {acc}. Either its gene (value found: {gene}) or product (Value found: {product}) was not found in the database genbank.\nPlease make sure this entry is from NCBI Protein db and it has the Gene/Protein information properly formated.")
    else:
        print(f">NCBI_PROTEIN~~~{gene}~~~{acc}~~~{product}~~~[{organism}]\n{seq}")

# close file connection
input_handle.close()
