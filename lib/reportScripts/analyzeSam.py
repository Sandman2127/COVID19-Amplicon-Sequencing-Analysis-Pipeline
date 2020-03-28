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
parser.add_argument('-threshold',type=str,default='5,5000,10,30',help='The default threshold to consider COVID + status. CSV format: <total genome sites required (int)>,<minimum reads at each amplicon site (int # of reads)>,<acceptable deviation from start site alignment (in bp)>,<MAQ required (default:30)>')
parser.add_argument('-bed',type=str,help='Expected alignment positions of the primers in the COVID19 genome in bed format (starting coordinate of the read is all that is necessary)')

args = parser.parse_args()
SAMF=args.inputSam
PATH=os.path.dirname(SAMF)
if PATH == "":
    PATH = "."
basename=os.path.basename(SAMF)[:-4] 
threshold=args.threshold
bedFile=args.bed

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
    plt.legend(loc='lower left', labelspacing=0.5, bbox_to_anchor= (1.04, 0.5), borderaxespad=0, frameon=False)
    plt.tick_params(axis='x', which='major', labelsize=8)
    plt.ylabel("Reads Aligned @ Map Quality Score")
    #plt.tight_layout()
    output_PNG = PATH + "/Amplicon_Barplot.png"
    plt.savefig(output_PNG,format="png",figsize=(12,8),bbox_inches='tight', dpi=300)

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
    output_PNG = PATH + "/Amplicon_Pos_Hist.png"
    fig.savefig(output_PNG,format="png", dpi=300)

def predCovid19Status(threshold,bedFile):
    BEDDF=pd.read_csv(bedFile,sep="\t",header=None,index_col=False,names=['AMPLICON_ALIGNING_CONTIG','START','END'])
    # test.bed
    # lcl|NC_045512.2_gene_1  3921    4071
    # lcl|NC_045512.2_gene_10 421     571
    # lcl|NC_045512.2_gene_3  281     431
    # lcl|NC_045512.2_gene_7  71      221
    # lcl|NC_045512.2_gene_9  71      221

    # split threshold string:
    totalSites,RequiredReads,DistanceThreshold,MAQ=threshold.split(",")
    print(BEDDF.head())
    COVID_prediction_requirements = "Your chosen parameters require: " + totalSites + " total amplicon sites be enriched for >= " +  RequiredReads + " high quality aligning reads at MAQ of >= " + MAQ + " with a distance threshold of +/- " + DistanceThreshold + " bp to consider a diagnosis + for COVID19."
    print(COVID_prediction_requirements)

    # perform filtering of reads for the above parameters:
    diagnostic_DFs = subSelectDF[subSelectDF.MAPQ >= int(MAQ)].groupby('RNAME')
    # set total accepted amplicons:
    amplicons_accepted = 0
    amplicon_list=[]
    for DF in range(0,len(diagnostic_DFs)):
        workingDF=list(diagnostic_DFs)[DF][1]
        workingDF.reset_index(drop=True)
        # pick out amplicon sequence name:
        amplicon=str(workingDF['RNAME'].iloc[0]).strip()
        # identify the amplicon expected start position:
        for index, row in BEDDF.iterrows():
            if str(row['AMPLICON_ALIGNING_CONTIG']) == amplicon:
                startPos = int(row['START'])
                break
        # set parameters for acceptable alignments:
        acceptedStartPosUpstream = startPos - int(DistanceThreshold)
        acceptedStartPosDownstream = startPos + int(DistanceThreshold)
        if acceptedStartPosUpstream < 0:
            acceptedStartPosUpstream = 0
        # test for aligning reads within expected position:
        Criterion_1 = workingDF.POS >= acceptedStartPosUpstream
        Criterion_2 = workingDF.POS <= acceptedStartPosDownstream
        AllCriterion = Criterion_1 & Criterion_2
        total_passing_reads = len(workingDF[AllCriterion].index)
        # append appropriate output:
        amplicon_data = amplicon + ":" + str(total_passing_reads)
        amplicon_list.append(amplicon_data)
        # increment total amplicons expected if total passing reads is greater than required_read threshold:
        if total_passing_reads > int(RequiredReads):
            amplicons_accepted = amplicons_accepted + 1
    
    # make diagnosis: 
    if amplicons_accepted >= int(totalSites):
        diagnosis = "+"
    else:
        diagnosis = "-"
    
    # write simple CSV output file
    fout = open("sample_summary.txt","w+")
    output = basename + "," + diagnosis + "," + "amplicons_accepted:" + str(amplicons_accepted) + "," + "amplicons_required:" + str(totalSites)
    for amp in amplicon_list:
        output = output + "," + amp
    fout.write(output)
    fout.close()

def main():
    print("Plotting alignment to all chromosomes x MAPQ...")
    plotBarByAmplicon()
    print("Plotting alignment start position for all reads on each chromosome...")
    plotHistByAmpliconPos()
    print("Predicting COVID status...")
    predCovid19Status(threshold,bedFile)
    print("Analysis Complete")

if __name__ == "__main__":
    print("Executing Main")
    main() 
