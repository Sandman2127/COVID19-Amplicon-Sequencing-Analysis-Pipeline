#!/usr/bin/env nextflow
 
/*
 * Defines the pipeline inputs parameters (giving a default value for each for them)
 * Each of the following parameters can be specified as command line options
 */

genome = file(params.genome)
barcodes = file(params.barcodes)
input_ch = Channel.fromPath(params.inputF)

lib = "$baseDir/lib"
reportScripts = "$lib/reportScripts"

process demuxSamples {
    input:
    file fastq from input_ch

    output:
    path '*.fq' into demux_ch

    //#Barcodes file
    //CTCTCCAG  RED_CLOUD_001.fq

    """
    sabre se -f $fastq -b $barcodes -u unknown_barcode.fastq > sabre.log.txt
    """
}


process alignBowtie2 {
    input:
    file demux_fq from demux_ch.flatMap()

    output:
    file "${demux_fq}.bam" into aligned_bam_sequences
    file "${demux_fq}.sam" into aligned_sam_sequences

    """
    # very senstive end to end alignment with 1 possible mismatch in the seed into a samtools view for only high quality alignments, followed by sorting and indexing
    bowtie2 --threads 4 --very-sensitive -N 1 -x $genome $demux_fq | samtools view -bh -q 5 - | samtools sort -@ 2 - > ${demux_fq}.bam ;

    # index alignments
    samtools index -@ 4 ${demux_fq}.bam ;

    # generate sam file for plotting analysis
    samtools view -h ${demux_fq}.bam > ${demux_fq}.sam
    """
}


process plotSamOutput {
    input:
    file sam from aligned_sam_sequences

    
    """
    python $reportScripts/analyzeSam.py -inputSam $sam 
    """

}