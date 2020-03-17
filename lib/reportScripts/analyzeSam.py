#!/usr/bin/env python3
  
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import numpy as np
import argparse
import os, sys
plt.switch_backend('agg')


parser = argparse.ArgumentParser(description='Generate plots of sam statistics')
parser.add_argument('-inputSam', type=str, help='The absolute path to the your input.sam, i.e. /home/Desktop/yourinput.sam', required=True)

args = parser.parse_args()
SAMF=args.inputSam
PATH=os.path.dirname(SAMF)
if PATH == "":
    PATH = "."
basename=os.path.basename(SAMF)[:-4] 

SAMDF = pd.read_csv(SAMF,sep='\s+',index_col=False,skiprows=13,names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL','VAL','VAL2','VAL3','VAL4','VAL5','VAL6','VAL7'],header=None)
#amplicon0-4166572      0       lcl|NC_045512.2_gene_1  3921    42      150M    *       0       0       TGGCAGTACTGAAATGCTAGCGAAATCTCTGAGAAAAGTGCCAACAGACAATTATATAACCACTTACCCGGGTCAGGGTTTAAATGGTTACACTGTAGAGGAGGCAAAGACAGTGCTTAAAAAGTGTAAAAGTGCCTTTTACATTCTACC  FG1G##GGB*G9G#G90GGGGGG4G#DG#GGGGGFG9GGGCGGGGG+GCGGGGGG8GFGG@G;GGGGGG*GGGGG/CGGFAF:GGGGFGFGGGGCDGGGGGGGGE<GGGGGGCGEGGGD+CGGGGGGG<GGGFGGGGGGGFG?G@FEG7G  AS:i:-6 XN:i:0  XM:i:3XO:i:0    XG:i:0  NM:i:3  MD:Z:5C19G2T121 YT:Z:UU
subSelectDF = SAMDF.loc[:,'QNAME':'TLEN']


def plotByAmplicon():
    MQgt40 = subSelectDF[subSelectDF.MAPQ >= 40 ].groupby('RNAME').count().loc[:,'MAPQ']
    MQgt30 = subSelectDF[subSelectDF.MAPQ >= 30 ].groupby('RNAME').count().loc[:,'MAPQ']
    MQgt20 = subSelectDF[subSelectDF.MAPQ >= 20 ].groupby('RNAME').count().loc[:,'MAPQ']
    MQgt10 = subSelectDF[subSelectDF.MAPQ >= 10 ].groupby('RNAME').count().loc[:,'MAPQ']
    merged = pd.merge(MQgt10, MQgt20, on='RNAME').merge(MQgt30,on='RNAME').merge(MQgt40,on='RNAME')
    merged.columns = ['MAQ>10','MAQ>20','MAQ>30','MAQ>40']
    merged.plot.bar()
    output_PNG = PATH + "/" + basename + "_Amplicon_Barplot.png"
    output_EPS = PATH + "/" + basename + "_Amplicon_Barplot.eps"
    plt.savefig(output_PNG,format="png", dpi=600)
    plt.savefig(output_EPS,format="eps", dpi=600)


def main():
    print("Plotting Alignment to all chromosomes x MAPQ")
    plotByAmplicon()


if __name__ == "__main__":
    print("Executing Main")
    main()