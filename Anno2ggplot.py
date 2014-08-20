#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse

def usage():
    test="name"
    message='''
python Anno2ggplot.py --input chr6.headingdate.annotation

Rank:
0	No variants
1	DOWNSTREAM;INTRON
2	SYNONYMOUS_CODING
3	UPSTREAM
4	NON_SYNONYMOUS_CODING

input
LOC_Os06g04960	DOWNSTREAM;INTRON	INTRON	NA	NA
output
LOC_Os06g04960	1	1	0	0

    '''
    print message

'''find unit of eff in rank, if found return 1 else return 0'''
def highimpact(eff, rank):
    flag = 0
    for e in eff:
        e = e.upper()
        if e in rank:
            flag = 1
            return flag
    return flag

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')    
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    high = ['SPLICE_SITE_ACCEPTOR',
'SPLICE_SITE_DONOR',
'START_LOST',
'EXON_DELETED',
'FRAME_SHIFT',
'STOP_GAINED',
'STOP_LOST',
'RARE_AMINO_ACID']
    moderate = ['NON_SYNONYMOUS_CODING',
'CODON_CHANGE',
'CODON_INSERTION',
'CODON_CHANGE_PLUS_CODON_INSERTION',
'CODON_DELETION',
'CODON_CHANGE_PLUS_CODON_DELETION',
'UTR_5_DELETED',
'UTR_3_DELETED',]
    low = ['SYNONYMOUS_START',
'NON_SYNONYMOUS_START',
'START_GAINED',
'SYNONYMOUS_CODING',
'SYNONYMOUS_STOP',]
    modifier = ['UTR_5_PRIME',
'UTR_3_PRIME',
'REGULATION',
'UPSTREAM',
'DOWNSTREAM',
'GENE',
'TRANSCRIPT',
'EXON',
'INTRON_CONSERVED',
'INTRON',
'INTRAGENIC',
'INTERGENIC',
'INTERGENIC_CONSERVED',
'NONE',
'CHROMOSOME',
'CUSTOM',
'CDS',]
    print 'Gene\tSNPeff\tINDELeff\tDELeff\tINTeff'
    with open (args.input, 'r') as infile:
        head = infile.readline()
        for line in infile:
            line = line.rstrip()
            if len(line) < 2:
                continue
            unit = re.split(r'\t',line)
            outline = []
            outline.append(unit[0])
            for i in range(1, 5):
                eff = re.split(r';', unit[i])
                if highimpact(eff, high):
                    outline.append('4') 
                elif highimpact(eff, moderate):
                    outline.append('3')
                elif highimpact(eff, low):
                    outline.append('2')
                elif highimpact(eff, modifier):
                    outline.append('1')
                else:
                    outline.append('0')
            print '\t'.join(outline)

if __name__ == '__main__':
    main()

