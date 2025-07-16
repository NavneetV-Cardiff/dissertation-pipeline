#!/usr/bin/env nextflow

// Define parameters
params.genome_url = "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.22.fa.gz"
params.fastq1 = "ERR166333_1.fastq"
params.fastq2 = "ERR166333_2.fastq"

// Define processes
process DOWNLOAD_GENOME {
    output:
    path "genome.fa.gz"

    script:
    """
    wget -O genome.fa.gz ${params.genome_url}
    gunzip genome.fa.gz
    """
}

process INDEX_GENOME {
    input:
    path "genome.fa"

    output:
    path "genome.fa.*"

    script:
    """
    bwa index -a bwtsw genome.fa
    """
}

process TRIM_FASTQ {
    input:
    path params.fastq1
    path params.fastq2

    output:
    path "trimmed_1.fastq.gz"
    path "trimmed_2.fastq.gz"

    script:
    """
    java -jar /apps/genomics/trimmomatic/0.39/trimmomatic-0.39.jar PE -phred33 \
        ${params.fastq1} ${params.fastq2} \
        trimmed_1.fastq.gz unpaired_1.fastq.gz \
        trimmed_2.fastq.gz unpaired_2.fastq.gz \
        ILLUMINACLIP:TruSeq-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// Define workflow
workflow {
    DOWNLOAD_GENOME()
    INDEX_GENOME()
    TRIM_FASTQ()
}