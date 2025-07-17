#!/usr/bin/env nextflow

// Define parameters
params.sra_id = "ERR166333"
params.genome_url = "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
params.chromosomes = [
    "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa.gz",
    "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.3.fa.gz",
    "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.4.fa.gz",
    "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.5.fa.gz",
    "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.9.fa.gz"
]
params.known_sites = "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
params.supporting_vcf = "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf"

// Define processes

// Step 1: Download SRA data and convert to FASTQ
process DownloadSRA {
    input:
    val params.sra_id

    output:
    path "${params.sra_id}_1.fastq.gz"
    path "${params.sra_id}_2.fastq.gz"

    script:
    """
    module load sra-tools/2.10.0
    prefetch ${params.sra_id}
    fasterq-dump ${params.sra_id}
    gzip ${params.sra_id}*.fastq
    """
}

// Step 2: Download and prepare the genome
process DownloadGenome {
    input:
    val params.genome_url
    val params.chromosomes

    output:
    path "Homo_sapiens.GRCh38.dna_sm.slim.fa"

    script:
    """
    wget -O genome.fa.gz ${params.genome_url}
    gunzip genome.fa.gz

    for url in ${params.chromosomes}; do
        wget \$url
    done

    gunzip *.gz
    cat *chromosome*fa > Homo_sapiens.GRCh38.dna_sm.slim.fa
    """
}

// Step 3: Index the genome
process IndexGenome {
    input:
    path "Homo_sapiens.GRCh38.dna_sm.slim.fa"

    output:
    path "Homo_sapiens.GRCh38.dna_sm.slim.fa.*"

    script:
    """
    module load bwa/0.7.17
    bwa index -a bwtsw Homo_sapiens.GRCh38.dna_sm.slim.fa
    """
}

// Step 4: Trim FASTQ files
process TrimFastq {
    input:
    path "${params.sra_id}_1.fastq.gz"
    path "${params.sra_id}_2.fastq.gz"

    output:
    path "${params.sra_id}.trim_1.fastq.gz"
    path "${params.sra_id}.trim_2.fastq.gz"

    script:
    """
    module load java/1.8
    module load trimmomatic/0.39
    java -jar /apps/genomics/trimmomatic/0.39/trimmomatic-0.39.jar PE -phred33 \
        ${params.sra_id}_1.fastq.gz ${params.sra_id}_2.fastq.gz \
        ${params.sra_id}.trim_1.fastq.gz unpaired_1.fastq.gz \
        ${params.sra_id}.trim_2.fastq.gz unpaired_2.fastq.gz \
        ILLUMINACLIP:TruSeq-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// Step 5: Map reads to the genome
process MapReads {
    input:
    path "Homo_sapiens.GRCh38.dna_sm.slim.fa"
    path "${params.sra_id}.trim_1.fastq.gz"
    path "${params.sra_id}.trim_2.fastq.gz"

    output:
    path "${params.sra_id}.sam"

    script:
    """
    module load bwa/0.7.17
    bwa mem Homo_sapiens.GRCh38.dna_sm.slim.fa ${params.sra_id}.trim_1.fastq.gz ${params.sra_id}.trim_2.fastq.gz > ${params.sra_id}.sam
    """
}

// Step 6: Convert SAM to BAM and sort
process SortBam {
    input:
    path "${params.sra_id}.sam"

    output:
    path "${params.sra_id}.sorted.bam"

    script:
    """
    module load samtools/1.17
    samtools view -bS ${params.sra_id}.sam | samtools sort -o ${params.sra_id}.sorted.bam
    """
}

// Step 7: Mark duplicates
process MarkDuplicates {
    input:
    path "${params.sra_id}.sorted.bam"

    output:
    path "${params.sra_id}.sorted.markdup.bam"

    script:
    """
    module load java/1.8
    module load picard/2.27.5
    java -jar \$PICARD MarkDuplicates I=${params.sra_id}.sorted.bam O=${params.sra_id}.sorted.markdup.bam M=${params.sra_id}.metrics REMOVE_DUPLICATES=false
    """
}

// Step 8: Base recalibration
process BaseRecalibration {
    input:
    path "${params.sra_id}.sorted.markdup.bam"
    path "Homo_sapiens.GRCh38.dna_sm.slim.fa"
    val params.known_sites

    output:
    path "${params.sra_id}.recal_data.table"

    script:
    """
    module load GATK/4.1.2.0
    gatk BaseRecalibrator -I ${params.sra_id}.sorted.markdup.bam -R Homo_sapiens.GRCh38.dna_sm.slim.fa --known-sites ${params.known_sites} -O ${params.sra_id}.recal_data.table
    """
}

// Step 9: Apply recalibration
process ApplyBQSR {
    input:
    path "${params.sra_id}.sorted.markdup.bam"
    path "${params.sra_id}.recal_data.table"
    path "Homo_sapiens.GRCh38.dna_sm.slim.fa"

    output:
    path "${params.sra_id}.recal.bam"

    script:
    """
    module load GATK/4.1.2.0
    gatk ApplyBQSR -R Homo_sapiens.GRCh38.dna_sm.slim.fa -I ${params.sra_id}.sorted.markdup.bam --bqsr-recal-file ${params.sra_id}.recal_data.table -O ${params.sra_id}.recal.bam
    """
}

// Step 10: Call variants
process CallVariants {
    input:
    path "${params.sra_id}.recal.bam"
    path "Homo_sapiens.GRCh38.dna_sm.slim.fa"

    output:
    path "${params.sra_id}.vcf.gz"

    script:
    """
    module load GATK/4.1.2.0
    gatk --java-options "-Xmx4g" HaplotypeCaller -R Homo_sapiens.GRCh38.dna_sm.slim.fa -I ${params.sra_id}.recal.bam -O ${params.sra_id}.vcf.gz
    """
}

// Step 11: Annotate variants
process AnnotateVariants {
    input:
    path "${params.sra_id}.vcf.gz"
    path "Homo_sapiens.GRCh38.dna_sm.slim.fa"

    output:
    path "${params.sra_id}.annot.vcf.gz"

    script:
    """
    module load GATK/4.1.2.0
    gatk Funcotator --variant ${params.sra_id}.vcf.gz --reference Homo_sapiens.GRCh38.dna_sm.slim.fa --ref-version hg38 --data-sources-path funcotator_dataSources.v1.6.20190124g --output ${params.sra_id}.annot.vcf.gz --output-file-format VCF
    """
}

// Define the workflow
workflow {
    DownloadSRA(params.sra_id)
    DownloadGenome(params.genome_url, params.chromosomes)
    IndexGenome()
    TrimFastq()
    MapReads()
    SortBam()
    MarkDuplicates()
    BaseRecalibration()
    ApplyBQSR()
    CallVariants()
    AnnotateVariants()
}