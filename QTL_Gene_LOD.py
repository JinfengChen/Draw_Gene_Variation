#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
import pandas as pd

def usage():
    test="name"
    message='''
python QTL_Gene_LOD.py --gff HeadingDate1_Chr3_962389_1665631.gff --lod ../input/MPR.cross.uniq.QTL.mr.table.new --trait HeadingDays --project HeadingDate1_Chr3_962389_1665631

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

'''
Chr6    MSU_osa1r7      gene    1274170 1274692 .       +       .       ID=LOC_Os06g03330;Name=LOC_Os06g03330;Note=hypothetical%20protein
Chr6    MSU_osa1r7      mRNA    1274170 1274692 .       +       .       ID=LOC_Os06g03330.1;Name=LOC_Os06g03330.1;Parent=LOC_Os06g03330
'''
def readgff(infile):
    data = defaultdict(str)
    s = re.compile(r'ID=(.*?);Name')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if line.startswith('Chr'): 
                unit = re.split(r'\t',line)
                if unit[2] == 'gene':
                    m = s.search(unit[8])
                    gene = m.groups(0)[0] if m else 'NA'
                    if not gene == 'NA':
                        data[gene] = [unit[0], unit[3], unit[4]]
    return data

def lodfile(lodfile, genelist, trait, prefix):
    lod_data = pd.read_table(lodfile)
    traits = re.split(r',', trait)
    data = defaultdict(lambda : str)
    for t in traits:
        temp = t.replace(r'/','')
        outfile = '%s.LOD' %(prefix)
        ofile = open(outfile, 'w')
        print >> ofile, 'Gene\tChr\tGene_start\tGene_end\tBin_start\tBin_end\tLOD'
        for g in sorted(genelist.keys()):
            chrs = genelist[g][0]
            start= int(genelist[g][1])
            end  = int(genelist[g][2])
            for i in range(len(lod_data['Position'])):
                #print '%s\t%s\t%s\t%s' %(g, chrs, start, end)
                #print '%s\t%s\t%s\t%s' %(str(i), lod_data['chr'][i], lod_data['Bin_start'][i], lod_data['Bin_end'][i])
                if lod_data['chr'][i] == chrs and int(lod_data['Bin_start'][i]) < start and int(lod_data['Bin_end'][i]) > end:
                    #print '%s\t%s\t%s\t%s' %(g, chrs, start, end)
                    #print '%s\t%s\t%s\t%s' %(str(i), lod_data['chr'][i], lod_data['Bin_start'][i], lod_data['Bin_end'][i])
                    line = '%s\t%s\t%s\t%s\t%s\t%s\t%s' %(g, chrs, start, end, lod_data['Bin_start'][i], lod_data['Bin_end'][i], lod_data[t][i])
                    data[g] = line
                elif lod_data['chr'][i] == chrs and int(lod_data['Bin_start'][i]) < start and int(lod_data['Bin_end'][i]) > start:
                    line = '%s\t%s\t%s\t%s\t%s\t%s\t%s' %(g, chrs, start, end, lod_data['Bin_start'][i], lod_data['Bin_end'][i], lod_data[t][i])
                    data[g] = line
                elif lod_data['chr'][i] == chrs and int(lod_data['Bin_start'][i]) < end and int(lod_data['Bin_end'][i]) > end:
                    line = '%s\t%s\t%s\t%s\t%s\t%s\t%s' %(g, chrs, start, end, lod_data['Bin_start'][i], lod_data['Bin_end'][i], lod_data[t][i])
                    data[g] = line
            if data.has_key(g):
                print >> ofile, data[g]
            else:
                print >> ofile, '%s\t%s\t%s\t%s\tNA\tNA\tNA' %(g, chrs, start, end)
        ofile.close()
            
def sublodfile(lodfile, chrs, start, end, qtl_trait, prefix):
    outfile = prefix + '.LOD.subfile' 
    ofile = open(outfile, 'w')
    print >> ofile, 'Position\tBin_start\tBin_end\tChr\tLOD'
    lod_data = pd.read_table(lodfile)
    position = []
    for i in range(len(lod_data['Position'])):
        #print lod_data['chr'][i], lod_data['Position'][i], chrs, start, end
        if lod_data['chr'][i] == chrs and int(lod_data['Position'][i]) > start and int(lod_data['Position'][i]) < end:
            line = '%s\t%s\t%s\t%s\t%s' %(lod_data['Position'][i], lod_data['Bin_start'][i], lod_data['Bin_end'][i], lod_data['chr'][i], lod_data[qtl_trait][i])
            position.append(lod_data['Bin_start'][i])
            position.append(lod_data['Bin_end'][i])
            print >> ofile, line
    ofile.close()
    position = sorted(position, key=int)
    start = position[0]
    end   = position[-1]
    return [outfile], start, end
 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gff')
    parser.add_argument('-l', '--lod')
    parser.add_argument('-t', '--trait')
    parser.add_argument('-p', '--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.gff) > 0
    except:
        usage()
        sys.exit(2)
    
    genelist = readgff(args.gff)
    lodfile(args.lod, genelist, args.trait, args.project)
    

if __name__ == '__main__':
    main()

