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


def plotBarByAmplicon():
    MQgt40 = subSelectDF[subSelectDF.MAPQ >= 40 ].groupby('RNAME').count().loc[:,'MAPQ']
    MQgt30 = subSelectDF[subSelectDF.MAPQ >= 30 ].groupby('RNAME').count().loc[:,'MAPQ']
    MQgt20 = subSelectDF[subSelectDF.MAPQ >= 20 ].groupby('RNAME').count().loc[:,'MAPQ']
    MQgt10 = subSelectDF[subSelectDF.MAPQ >= 10 ].groupby('RNAME').count().loc[:,'MAPQ']
    MQgt0 = subSelectDF.groupby('RNAME').count().loc[:,'MAPQ']
    merged = pd.merge(MQgt0, MQgt10, on='RNAME').merge(MQgt20,on='RNAME').merge(MQgt30,on='RNAME').merge(MQgt40,on='RNAME')
    merged.columns = ['MAQ>=0','MAQ>=10','MAQ>=20','MAQ>=30','MAQ>=40']
    merged.plot.bar()
    output_PNG = PATH + "/" + basename + "_Amplicon_Barplot.png"
    output_EPS = PATH + "/" + basename + "_Amplicon_Barplot.eps"
    plt.savefig(output_PNG,format="png", dpi=600)
    plt.savefig(output_EPS,format="eps", dpi=600)

def plotHistByAmpliconPos():
    #Group by each contig an amplicon aligns to
    groupByRNAME = subSelectDF.groupby('RNAME')
    lengthGBR=len(groupByRNAME)
    #plot figure
    fig, axs = plt.subplots(lengthGBR,sharex=False,constrained_layout=True,figsize=(8.5,12))
    fig.suptitle('Amplicon read alignment by contig and start position', fontsize=16)
    for num in range(0,lengthGBR):
        name =list(groupByRNAME)[num][0]
        df =list(groupByRNAME)[num][1]
        axs[num].hist(df.loc[:,'POS'])
        xlab = "Alignment position on contig: " + str(name)
        axs[num].set(xlabel=xlab,ylabel="Total Reads")
    output_PNG = PATH + "/" + basename + "_Amplicon_Pos_Hist.png"
    output_EPS = PATH + "/" + basename + "_Amplicon_Pos_Hist.eps"
    fig.savefig(output_PNG,format="png", dpi=600)
    fig.savefig(output_EPS,format="eps", dpi=600)


def main():
    print("Plotting alignment to all chromosomes x MAPQ")
    plotBarByAmplicon()
    print("Plotting alignment start position for all reads on each chromosome")
    plotHistByAmpliconPos()


if __name__ == "__main__":
    print("Executing Main")
    main()