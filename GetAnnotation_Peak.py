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
python GetAnnotation.py --input QTL.regions.list
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

'''
Chr	Start	End	Name	Length
Chr3	962389	1479666	HeadingDate1	517278
Chr6	1989540	2276708	HeadingDate3	287169
Chr8	26497269	26759802	HeadingDate3	262534
'''
def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if unit[1].isdigit():
                    print line
                    #awk '$1 ~ /Chr6/ && $4 > 2458650 && $5 < 2664019' ../input/MSU7.gene.2.gff > chr6.Rankcold.gff
                    prefix = '%s_%s_%s_%s' %(unit[3], unit[0], unit[1], unit[2])
                    cmd1 = 'awk \'$1 ~ /%s/ && $4 > %s && $5 < %s\' ../input/MSU7.gene.2.gff > %s.gff' %(unit[0], unit[1], unit[2], prefix)
                    cmd2 = 'perl AnnoGene.pl --gff %s.gff > %s.annotation' %(prefix, prefix)
                    cmd3 = 'python Anno2ggplot.py --input %s.annotation > %s.table' %(prefix, prefix)
                    cmd4 = 'cut -f1,12 %s.annotation > %s.annotation.brief' %(prefix, prefix)
                    cmd5 = 'python QTL_Gene_RNAseq.py --input %s.table --RNAseq ../input/QTL.ColdStress.FPKM.table' %(prefix)
                    cmd6 = 'python QTL_Gene_RNAseq.py --input %s.table --RNAseq ../input/QTL.Drought.FPKM.table' %(prefix)
                    cmd7 = 'python QTL_Gene_RNAseq.py --input %s.table --RNAseq ../input/QTL.Salt.FPKM.table' %(prefix)
                    cmd8 = 'python QTL_Peak_Gene_heatmap.py --input %s.table --rnaseq1 %s.ColdStress.table --rnaseq2 %s.Drought.table --rnaseq3 %s.Salt.table --anno %s.annotation.brief --output %s' %(prefix, prefix, prefix, prefix, prefix, prefix)
                    os.system(cmd1)
                    os.system(cmd2)
                    os.system(cmd3)
                    os.system(cmd4)
                    os.system(cmd5)
                    os.system(cmd6)
                    os.system(cmd7)
                    os.system(cmd8)
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    readtable(args.input)

if __name__ == '__main__':
    main()

