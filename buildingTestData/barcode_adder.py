#!/usr/bin/env python
#use python 3

import argparse

parser = argparse.ArgumentParser(description='Append Barcodes to Samfile')
parser.add_argument('-sam', type=str, help='sam from amplicon sequencing simulation ART', required=False)

args = parser.parse_args()
samF = args.sam

#amplicon0-4166666	0	amplicon0	1	99	2=1X11=1X7=1X3=1X7=1X115=	*	0	0	TGACACTACTGAAAAGCTAGCGCAAGGTTTGAGACAAGTGCCAACAGACAATTATATAACCACTTACCCGGGTCAGGGTTTAAATGGTTACACTGTAGAGGAGGCAAAGACAGTGCTTAAAAAGTGTAAAAGTGCCTTTTACATTCTACC	G5#GG7*#GGGGGG#GGGFG#E#CGF#GG#GGAG#1GE8G:GCGGFGC@GGCGGGEGG8EGGGGCGG*CGGFGGD>GGFAGCFCCGGGGGGGGGGDGGDGGGGGF7GEGGFGGGGGGGGGGFGGGGGFCGGGGGCGGGGGEGGGGGGGGG

barcodes = ['CTCTCCAG','TAATTG','ATCTCGT','GACAACT','CTCGCAA','TGGACACT','TGTCAAT','TCCTGCT','GAACTT','ATGCT','ATTCCAA','GACACACT','CGCGT','CATACGCG','CTATCACT','CTGAACCA','TCTCCGT','TGTACA','AAGCAACT','ACCGA','GTAAG','TGATCGCT','TGCGG','ACTAA','GAGGTCCT','TAGCTAT','CAGCGCAAGA','GCTCGCCAT','TGTACCAG','TGTACGCA','TTGGCGCT','GTTCACA','CATGG','ACTACAAT','GACTAACT','ATGGTGA','TATTGCAG','ATCTGACT','GTCACGA','AACGACCACA','CGCCTCAT','CTTATG','TAGAG','GGCAT','CCGACG','TGGTCAAG','ACCAAG','CCATCCAA','GTTCGGT','GCCGCAAT','CATAAG','TTGAGACAG','ACCGTCCAT','GCGTGCCAGA','CCGAT','TCCTCCA','ACACG','CGCAAGA','ACACAACA','ATATT','GTCTCAACG','CCGCA','TCGTGACAGT','AATTG','TCCGT','TATAAGCAG','ATTCA','ACATGCCAG','TGCCTA','AAGGCCAACT','ACTCCACG','GGTTG','TTCTCA','CTGCCGT','TTCCA','GAGCGCT','TAATTAA','TGTGAGG','TGTTGACG','TACCT','CCAGGA','GGATGA','ACAGAAT','ATACTGAG','CTCCAA','TTAGGA','CCAAGACAGT','CATTGA','TCATT','GAATAGA','TTCTG','ACCTAA','GCGTAG','CGTAGCAACA','AAGCAGA','CAATTGCT']

highQValue="G"
with open(samF,"r+") as FIN:
    lines = FIN.readlines()
    lineNo=0
    for line in lines:
        line = line.strip()
        if line[0] == "@":
            print(line)
        elif lineNo < 95:
            line = line.split()
            amp1,val,ampname,pos1,pos2,CIGAR,star,zero1,zero2,read,qualityScores = line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10]
            read = barcodes[lineNo] + str(read)
            lengthoBC = len(barcodes[lineNo])
            appendQuality = lengthoBC * highQValue
            qualityScores =  str(appendQuality) + qualityScores
            CIGAR = str(lengthoBC) + str("=") + str(CIGAR)
            outputstring = str(amp1) + "\t" + str(val) + "\t" +  str(ampname) + "\t" +  str(pos1) + "\t" +  str(pos2) + "\t" +  str(CIGAR) + "\t" + str(star) + "\t" + str(zero1) + "\t" +  str(zero2) + "\t" + str(read) + "\t" + str(qualityScores) 
            print(outputstring)
            lineNo = lineNo + 1
        else: # for line 95
            line = line.split()
            amp1,val,ampname,pos1,pos2,CIGAR,star,zero1,zero2,read,qualityScores = line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10]
            read = barcodes[lineNo] + str(read)
            lengthoBC = len(barcodes[lineNo])
            appendQuality = lengthoBC * highQValue
            qualityScores =  str(appendQuality) + qualityScores
            CIGAR = str(lengthoBC) + str("=") + str(CIGAR)
            outputstring = str(amp1) + "\t" + str(val) + "\t" +  str(ampname) + "\t" +  str(pos1) + "\t" +  str(pos2) + "\t" +  str(CIGAR) + "\t" + str(star) + "\t" + str(zero1) + "\t" +  str(zero2) + "\t" + str(read) + "\t" + str(qualityScores)
            print(outputstring)
            lineNo = 0
        





