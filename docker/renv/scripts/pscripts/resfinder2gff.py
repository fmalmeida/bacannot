#!/usr/bin/env python3
# coding: utf-8

## Def help message
"""
A simple script to be used as part of Bacannot pipeline. Created in order to aggregate the
Resfinder results to the main GFF output (which merges all annotations).
---
Copyright (C) 2020 Felipe Marques de Almeida (almeidafmarques@gmail.com)
License: Public Domain

Usage:
    resfinder2gff.py
    resfinder2gff.py -h|--help
    resfinder2gff.py -v|--version
    resfinder2gff.py [ --input <results_tab.txt> ]

Options:
    -h,--help                                    Show this screen.
    -v,--version                                 Show version information
    -i,--input=<results_tab.txt>                 Resfinder 'results_tab.txt' file
"""

##################################
### Loading Necessary Packages ###
##################################
from docopt import docopt
import pandas as pd
import re

####################################
### Def resfinder results loader ###
####################################
def load_resfinder(input):

    # Get results
    res_df = pd.read_csv(input, sep='\t', comment='#')
    res_df = res_df.sort_values(by=['Contig', 'Position in contig']).reset_index(drop=True)

    # Subset information
    for index, line in res_df.iterrows():

        num=f"Resfinder_{index+1}"
        gene=line['Resistance gene']
        id=line['Identity']
        contig=line['Contig']
        start=line['Position in contig'].split('..')[0]
        end=line['Position in contig'].split('..')[1]
        if start < end:
            strand="+"
        else:
            strand="-"
        source='Resfinder'
        type='Resistance'
        target=line['Phenotype'].replace(";","")
        acc=line['Accession no.']

        print(f"{contig}\t{source}\t{type}\t{start}\t{end}\t.\t{strand}\t.\tID={num};Additional_database={source};Resfinder_gene={gene};Resfinder_phenotype={target};Resfinder_reference={acc}".replace(" ", "_"))



############
### Main ###
############

if __name__ == '__main__':
    arguments = docopt(__doc__, version='v1.0 by Felipe Marques de Almeida')

    # Main pipe: resfinder to gff
    if arguments['--input']:
        load_resfinder(input=arguments['--input'])
