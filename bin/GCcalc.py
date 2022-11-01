#!/usr/bin/env python3

## calculate GC
from __future__ import division

# __author__ = "Wenchao Lin"
# __copyright__ = "Copyright 2014,Tianjin Biochip Corp."
# __license__ = "GPL"
# __version__ = "1.0.0"
# __maintainer__ = "Wenchao Lin"
# __email__ = "linwenchao@yeah.net"


from Bio import SeqIO
import sys
from optparse import OptionParser
import numpy as np

def GC_content_window(s):

    """
    Return GC content of input sequence
    """

    gc = sum(s.count(x) for x in ['G','C','g','c','S','s'])
    gc_content = gc/float(len(s))
    return round(gc_content,4)


def GC_skew_window(s):

    """
    Reuturn GC skew of a sequence
    """

    g = s.count('G')+s.count('g')
    c = s.count('C')+s.count('c')

    try:
        skew = (g-c)/float(g+c)
    except ZeroDivisionError:
        skew = 0
    return round(skew,4)


def main():

    """
    Calculate GC content and GC skew of input Fasta sequence
    """

    usage = "usage: %prog -f input.fa [-w 1000] [-s 1000]"
    parser = OptionParser(usage = usage)
    parser.add_option("-f","--file",dest="filename",help="Input Fasta format file",metavar="FASTA")
    parser.add_option("-w","--window",dest="WindowSize",help="default:1000 WindowSize to calculate",default=1000,type='int')
    parser.add_option("-s","--step",dest="StepSize",help="default:1000 StepSize for slide widows",default=1000,type='int')
    (options,args) = parser.parse_args()

    window = options.WindowSize
    step = options.StepSize
    seqobj = SeqIO.parse(options.filename,'fasta')

    for record in seqobj:
        pos_array = []
        gc_content_value_array = []
        gc_skew_value_array = []
        name = record.id
        seq = record.seq
        start = 0
        end = 0
        gc = 0
        gc_skew = 0

        for i in range(0,len(seq),step):
            subseq = seq[i:i+window]
            gc_content = (GC_content_window(subseq))
            gc_skew = (GC_skew_window(subseq))
            start = (i + 1 if (i+1<=len(seq)) else i)
            end = ( i + step if (i+ step<=len(seq)) else len(seq))
            print ("%s\t%s\t%s\t%s\t%s" % (name,start,end,gc_content,gc_skew))


if __name__ == '__main__':
    main()
    sys.exit()