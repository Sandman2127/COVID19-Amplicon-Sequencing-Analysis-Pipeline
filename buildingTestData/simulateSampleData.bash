#!/bin/bash

/home/miniconda/bin/art_illumina -ss MSv3 --samout -amp -na -i ./amplicons.fa -l 150 --fcov 4166666 -o single_ended_COVID_amplicon --rndSeed 127

python barcode_adder.py -sam single_ended_COVID_amplicon.sam > SE_result_barcodes.sam

singularity run /home/ec2-user/SIF/COVID19_Analysis.sif samtools fastq SE_result_barcodes.sam > SE_result_barcodes.fastq
