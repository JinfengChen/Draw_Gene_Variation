#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
python QTL_Gene_RNAseq.py --input chr3.headingdate.ggplot.table --RNAseq QTL.ColdStress.FPKM.table
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

'''
Gene	SNPeff	INDELeff	DELeff	INTeff
LOC_Os03g02762	0	1	0	0
LOC_Os03g02770	0	0	0	0
LOC_Os03g02780	0	0	0	0
LOC_Os03g02790	0	0	0	0
LOC_Os03g02800	0	0	0	0
LOC_Os03g02820	0	0	0	0
'''
def readgene(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                data[unit[0]] = line
    return data

def sub_rnaseq(gene, rnaseq, filename, stress):
    ofile = open(filename, 'w')
    print >> ofile, rnaseq['Gene']
    for g in sorted(gene.keys()):
        if g.startswith('LOC'):
            if rnaseq.has_key(g):
                print >> ofile, rnaseq[g]
            else:
                if stress == 'Salt':
                    print >> ofile, '%s\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0' %(g)
                else:
                    print >> ofile, '%s\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0' %(g) 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-R', '--RNAseq')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0 and len(args.RNAseq) > 0
    except:
        usage()
        sys.exit(2)
    
    gene   = readgene(args.input)
    rnaseq = readgene(args.RNAseq)

    #chr3.headingdate.ggplot.table
    prefix = os.path.splitext(str(args.input))[0] 
    #QTL.ColdStress.FPKM.table
    stress = re.split(r'\.', os.path.basename(str(args.RNAseq)))[1]

    filename = '%s.%s.table' %(prefix, stress)
    sub_rnaseq(gene, rnaseq, filename, stress)

if __name__ == '__main__':
    main()

