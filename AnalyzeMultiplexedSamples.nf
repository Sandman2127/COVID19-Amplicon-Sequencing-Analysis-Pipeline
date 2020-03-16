genome = file(params.genome)
barcodes = file(params.barcodes)
input_ch = Channel.fromPath(params.inputF)

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
    file "${demux_fq}.bam" into aligned_sequences_ch

    """
    bowtie2 --threads 4 --very-sensitive -N 1 -x $genome $demux_fq | samtools view -bh -q 30 - | samtools sort -@ 2 - > ${demux_fq}.bam ;
    # very senstive end to end alignment with 1 possible mismatch in the seed into a samtools view for only high quality alignments, followed by sorting and indexing

    samtools index -@ 4 ${demux_fq}.bam ;
    # index alignments
    """
}