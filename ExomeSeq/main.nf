#!/usr/bin/env nextflow

//this is the final version
params.fastq_id = ["ERR166333", "ERR166337"]
params.adapter_file = "/shared/home1/c.c24053373/pipeline/data/TruSeq-PE.fa"
params.out_dir = "/scratch/SCWF00079/shared/navneet/pipeline/output1"

//include { download_reference } from './modules/download_reference.nf'
//include { edit_ref_genome } from './modules/edit_ref_genome.nf'


process referenceDownload {
    tag "Downloading the reference genome"

    output: 
    path 'Homo_sapiens_assembly38.fasta', emit: ref_genome

    script:
    """
    wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
    """
}

process editRefGenome {

    container 'oras://community.wave.seqera.io/library/bwa_samtools_trimmomatic:2827ad58b1521516'
    input:
    path ref_genome

    output:
    path "combined_ref.fa", emit: chr_edited

    script:
    """
    samtools faidx ${ref_genome} \
    chr1 > combined_ref.fa
    """
}

process indexFastaBWA {

    container 'oras://community.wave.seqera.io/library/bwa_samtools_trimmomatic:2827ad58b1521516'

    input:
    path ref_genome

    output:
    tuple path(ref_genome), path ("${ref_genome}.*"), emit: reference_group

    script:
    """
    bwa index -a bwtsw $ref_genome
    """
}

process fastqTrim{

    tag "$fastq_id"

    container 'oras://community.wave.seqera.io/library/bwa_samtools_trimmomatic:2827ad58b1521516'

    
    input:
    tuple val (fastq_id), path(fastq_1), path (fastq_2)  // Paired-end FASTQ input files
    path adapter_file 

    output:
    tuple val(fastq_id), path("${fastq_id}.R1.trimmed.fastq.gz"), path("${fastq_id}.R2.trimmed.fastq.gz")

    script:
    """
    trimmomatic PE -phred33 \
    ${fastq_1} ${fastq_2} \
    ${fastq_id}.R1.trimmed.fastq.gz ${fastq_id}.R1.unpaired.fastq.gz \
    ${fastq_id}.R2.trimmed.fastq.gz ${fastq_id}.R2.unpaired.fastq.gz \
    ILLUMINACLIP:${adapter_file}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process performFastqc {

    tag "Performing fastqc on $fastq_id"
    publishDir params.out_dir, mode: 'copy', overwrite: true
    container 'oras://community.wave.seqera.io/library/fastqc:0.12.1--104d26ddd9519960'

    input:
    tuple val (fastq_id), path(fastq_1), path (fastq_2)

    output:
    path("*_fastqc.*")

    script:
    """
    fastqc ${fastq_1} ${fastq_2}
    """
}

process refGeneMapping {
    tag "mapping $fastq_id to reference genome"

    container 'oras://community.wave.seqera.io/library/bwa_samtools_trimmomatic:2827ad58b1521516'

    input:
    tuple val (fastq_id), path(fastq_1), path (fastq_2), path(ref_genome), path(reference_files)

    output:
    tuple val (fastq_id), path("${fastq_id}.sam")

    script:
    """
    bwa mem -R "@RG\\tID:${fastq_id}\\tSM:${fastq_id}\\tPL:ILLUMINA\\tLB:lib1" \
    ${ref_genome} ${fastq_1} ${fastq_2} > ${fastq_id}.sam

    """
}

process bamCreation {
    tag "making bam file for $fastq_id"
    container 'oras://community.wave.seqera.io/library/bwa_samtools_trimmomatic:2827ad58b1521516'

    input:
    tuple val (fastq_id), path(sam_file)

    output:
    tuple val (fastq_id), path("${fastq_id}.sorted.bam")

    script:
    """
    samtools view -bS ${sam_file} > ${fastq_id}.bam
    samtools sort -o ${fastq_id}.sorted.bam ${fastq_id}.bam 
    """
}

process markDuplicates {
    tag "Marking duplicates for $fastq_id"

    container 'oras://community.wave.seqera.io/library/picard:3.4.0--2976616e7cbd4840'

    input:
    tuple val (fastq_id), path (bam_file)

    output:
    tuple val (fastq_id), path("${fastq_id}.sorted.markdup.bam"), path("${fastq_id}.sorted.metrics")

    script:
    """
    picard MarkDuplicates \
    I=${bam_file} O=${fastq_id}.sorted.markdup.bam \
    M=${fastq_id}.sorted.metrics 
    REMOVE_DUPLICATES=false
    """

}

process bamStats {
    tag "$fastq_id"

    container 'oras://community.wave.seqera.io/library/bamtools:2.5.3--ad2f385a555ee16d'
    input: 
    tuple val(fastq_id), path (markdup_file), path(_)

    output:
    tuple val(fastq_id), path("${fastq_id}.sorted.markdup.stats")

    script:
    """
    bamtools stats -in $markdup_file > ${fastq_id}.sorted.markdup.stats
    """
}

process indexingBamFile {
    tag "performing indexingBamFile on $fastq_id"

    container 'oras://community.wave.seqera.io/library/bwa_samtools_trimmomatic:2827ad58b1521516'

    input: 
    tuple val(fastq_id), path (markdup_file), path(_)

    output:
    tuple val(fastq_id), path("${fastq_id}.sorted.markdup.bam.bai")

    script:
    """
    samtools index ${markdup_file}
    """
}

process addReadGroups {
    tag "Adding read groups to $fastq_id"

    container 'oras://community.wave.seqera.io/library/picard:3.4.0--2976616e7cbd4840'

    input:
    tuple val (fastq_id), path (markdup_file), path(_)

    output:
    tuple val(fastq_id), path("${fastq_id}.sorted.markup.rg.bam")

    script:
    """
    picard AddOrReplaceReadGroups \
    I=${markdup_file} O=${fastq_id}.sorted.markup.rg.bam SO=coordinate RGID=1 RGLB=libl RGPL=illumina RGPU=unit1 RGSM=${fastq_id} CREATE_INDEX=True
    """
}

process indexReference {
    tag "Indexing reference genome"
    container 'oras://community.wave.seqera.io/library/gatk_samtools:0465b74d4eaac4e4'

    input: 
    path ref_genome

    output:
    path "*.fai", emit: indexed_ref_fai

    script:

    """
    samtools faidx ${ref_genome}
    """
}

process downloadMillsand1000G {
    tag "Downloading Mills and 1000G gold standard indels VCF"
    
    output:
    path 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz', emit: var_reference_gz

    script:
    """
    wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz 
    """

}

process indexVariantFile {
    tag "Index variant reference panel"
    container '/shared/home1/c.c24053373/pipeline/singularity_cache/gatk_4.5.0.0.sif'

    input:
    path var_reference_gz

    output:
    tuple path(var_reference_gz), path('*.tbi'), emit: indexed_var_reference

    script:
    """
    gatk IndexFeatureFile -I ${var_reference_gz}
    """
}

process createSequenceDict {
    tag "Create sequence dictionary"

    container 'oras://community.wave.seqera.io/library/picard:3.4.0--2976616e7cbd4840'

    input:
    path ref_genome

    output:
    path 'combined_ref.dict', emit: reference_seq_dict

    script:
    """
    picard CreateSequenceDictionary \
    R=${ref_genome} \
    O=combined_ref.dict
    """
}

process createRecalModel {
    tag "$fasta_id"

    container '/shared/home1/c.c24053373/pipeline/singularity_cache/gatk_4.5.0.0.sif'

    input:
    tuple val(fasta_id), path(bam_with_readgroups), path(ref_genome), path(var_reference_gz), path(seq_dict), path(indexed_ref_fai), path(indexed_var_reference)

    output:
    tuple val(fasta_id), path("${fasta_id}.recal_data.table"), path(bam_with_readgroups), path(ref_genome), path(var_reference_gz), path(seq_dict), path(indexed_ref_fai), path(indexed_var_reference)

    script:
    """
    gatk BaseRecalibrator \
        -I $bam_with_readgroups \
        -R $ref_genome \
        --known-sites $var_reference_gz \
        -O ${fasta_id}.recal_data.table
    """
}

process recalBam {
    tag "$fastq_id"

    container '/shared/home1/c.c24053373/pipeline/singularity_cache/gatk_4.5.0.0.sif'

    input:
    tuple val(fastq_id), path(recal_file), path(bam_with_readgroups), path(ref_genome), path(indexed_var_reference), path(var_reference_gz), path(indexed_ref_fai), path(seq_dict)

    output:
    tuple val(fastq_id), path("${fastq_id}.sorted.markdup.rg.recal.bam"), path(ref_genome), path(indexed_ref_fai), path(seq_dict)

    script:
    """
    gatk ApplyBQSR -R $ref_genome -I $bam_with_readgroups --bqsr-recal-file $recal_file -O ${fastq_id}.sorted.markdup.rg.recal.bam
    """

}


process callVariants {
    tag "$fastq_id"

    container '/shared/home1/c.c24053373/pipeline/singularity_cache/gatk_4.5.0.0.sif'

    input:
    tuple val(fastq_id), path(recal_bam), path(ref_genome), path(indexed_ref_fai), path(seq_dict)

    output:
    tuple val(fastq_id), path("*.vcf.gz"), path(recal_bam), path(ref_genome), path(indexed_ref_fai), path(seq_dict)

    script:
    """
    gatk --java-options "-Xmx4g" HaplotypeCaller -R ${ref_genome} -I ${recal_bam} -O ${fastq_id}.vcf.gz    
    """
}


process downloadVariantArrays {
    tag "Download 1000 genome vcf"
    input:    

    output:
    path '1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf', emit: variant_array_vcf
    path '1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf.idx', emit: variant_array_idx

    script:
    """
    wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf
    wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf.idx
    """
}

process indexCalledVariants {
    tag "$fastq_id"

    container '/shared/home1/c.c24053373/pipeline/singularity_cache/gatk_4.5.0.0.sif'

    input:
    tuple val(fastq_id), path(vcf_file), path(recal_bam), path(ref_genome), path(indexed_ref_fai), path(seq_dict)

    output:
    tuple val(fastq_id), path(vcf_file), path('*.tbi'), emit: indexed_vcf

    script:
    """
    gatk IndexFeatureFile -I ${vcf_file}
    """
}

process purgeIndels {
    tag "$fastq_id"

    container '/shared/home1/c.c24053373/pipeline/singularity_cache/gatk_4.5.0.0.sif'

    input:
    tuple val(fastq_id), path(vcf_file), path(indexed_vcf)

    output:
    tuple val(fastq_id), path("*snp.vcf.gz"), path(vcf_file), path(indexed_vcf)

    script:
    """
    gatk SelectVariants --variant ${vcf_file} --select-type SNP --output ${fastq_id}.snp.vcf.gz
    """
}

process indexSNPS {
    tag "$fastq_id"

    container '/shared/home1/c.c24053373/pipeline/singularity_cache/gatk_4.5.0.0.sif'

    input:
    tuple val(fastq_id), path(snp_vcf), path(vcf_file), path(indexed_vcf)

    output:
    tuple val(fastq_id), path(snp_vcf), path('*.tbi'), emit: indexed_snp_vcf

    script:
    """
    gatk IndexFeatureFile -I ${snp_vcf}
    """
}

process filterVariants {
    tag "$fastq_id"

    container '/shared/home1/c.c24053373/pipeline/singularity_cache/gatk_4.5.0.0.sif'

    input:
    tuple val(fastq_id), path(snp_vcf), path(indexed_snp_vcf)

    output:
    tuple val(fastq_id), path("*.snp.filtered.vcf.gz"), path(snp_vcf), path(indexed_snp_vcf)

    script:
    """
    gatk VariantFiltration  --variant $snp_vcf --filter-expression "QD<2.0" --filter-name "QD2" --filter-expression "QUAL<30.0" --filter-name "QUAL30"  --filter-expression "SOR>3.0" --filter-name "SOR3" --filter-expression "FS>60.0" --filter-name "FS60" --filter-expression "MQ<40.0" --filter-name "MQ40" --filter-expression "MQRankSum<-12.5" --filter-name "MQRankSum-12.5" --filter-expression "ReadPosRankSum<-8.0" --filter-name "ReadPosRankSum-8" --output ${fastq_id}.snp.filtered.vcf.gz
    """
}

process indexVCF {
    tag "$fastq_id"

    container '/shared/home1/c.c24053373/pipeline/singularity_cache/gatk_4.5.0.0.sif'

    input:
    tuple val(fastq_id), path(filtered_vcf_file), path(vcf_file), path(indexed_vcf)

    output:
    tuple val(fastq_id), path(filtered_vcf_file), path('*.tbi'), emit: filtered_indexed_vcf

    script:
    """
    gatk IndexFeatureFile -I ${filtered_vcf_file}
    """
}

process calcGenotypePosteriors {
    tag "$fastq_id"

    container '/shared/home1/c.c24053373/pipeline/singularity_cache/gatk_4.5.0.0.sif'

    input:
    tuple val(fastq_id), path(snp_vcf), path(indexed_snp_vcf), path(thousand_genome_vcf), path(thousand_genome_vcf_idx)

    output:
    tuple val(fastq_id), path("*.snp.filtered.ref.vcf.gz"), path(thousand_genome_vcf), path(thousand_genome_vcf_idx), emit: refined_vcf

    script:
    """
    gatk --java-options "-Xmx4g" CalculateGenotypePosteriors \
    -V ${snp_vcf} \
    -O ${fastq_id}.snp.filtered.ref.vcf.gz \
    -supporting ${thousand_genome_vcf}
    """
}

process variantAnnotation {
    tag "$fastq_id"

    input:
    tuple val(fastq_id), path(refined_vcf_file), path(thousand_genome_vcf), path(thousand_genome_vcf_idx)

    output:
    path "*.hg38_multianno.txt"
    path "*.hg38_multianno.vcf"

    container '/scratch/SCWF00079/shared/navneet/pipeline/annovar_2018Apr16.sif'

    publishDir params.out_dir, mode: 'copy', overwrite: true

    script:
    """

    # Download necessary ANNOVAR databases (e.g., refGene, avsnp150, clinvar)
    #[ ! -f humandb/hg38_refGene.txt ] && annotate_variation.pl -downdb -buildver hg38 -webfrom annovar refGene humandb/
    #[ ! -f humandb/hg38_cytoBand.txt ] && annotate_variation.pl -build hg38 -downdb cytoBand humandb/
    #[ ! -f humandb/hg38_avsnp150.txt ] && annotate_variation.pl -downdb -buildver hg38 -webfrom annovar avsnp150 humandb/
    #[ ! -f humandb/hg38_clinvar.txt ] && annotate_variation.pl -downdb -buildver hg38 -webfrom annovar clinvar_20250721 humandb/

    # Convert VCF to ANNOVAR input format
    convert2annovar.pl -format vcf4 ${refined_vcf_file} > ${refined_vcf_file.baseName}.avinput

    # Annotate variants using ANNOVAR
        table_annovar.pl \
        ${refined_vcf_file} \
        /scratch/SCWF00079/shared/navneet/pipeline/humandb \
        -buildver hg38 \
        -out ${refined_vcf_file.baseName} \
        -remove \
        -protocol refGene,cytoBand,avsnp150,clinvar_20250721 \
        -operation g,r,f,f \
        -nastring . \
        -vcfinput
        

    """
}


workflow {

    ref_genome = referenceDownload()

    chr_edited = editRefGenome(ref_genome)

    indexed_reference = indexFastaBWA(chr_edited)

    adapter_ch = file(params.adapter_file)

    fastq_ch = Channel.fromFilePairs("data/*_{1,2}*.fastq.gz")
                      .map { sample_id, reads -> tuple(sample_id, reads[0], reads[1]) }

    fastq_ch.view { "FASTQ ID: ${it[0]}, Read 1: ${it[1]}, Read 2: ${it[2]}" }
    trim_output = fastqTrim(fastq_ch, adapter_ch)

    performFastqc(trim_output)

    mapping_input = trim_output.combine(indexed_reference)

    mapped_data = refGeneMapping(mapping_input)

    bam_file = bamCreation(mapped_data)

    markedup_files = markDuplicates(bam_file)

    bamStats_report = bamStats(markedup_files)

    indexed_bam = indexingBamFile(markedup_files)

    bam_with_readgroups = addReadGroups(markedup_files)

    reference_fai_file = indexReference(chr_edited)

    var_reference_gz = downloadMillsand1000G()

    indexed_var_reference = indexVariantFile(var_reference_gz)

    seq_dict = createSequenceDict(chr_edited)

    recalibration_inputs = bam_with_readgroups
        .combine(chr_edited.collect())        // add FASTA
        .combine(indexed_var_reference.collect())   // add known-sites VCF
        .combine(reference_fai_file.collect()) // add .fai
        .combine(seq_dict.collect())
        .map { fastq_id, bam, fasta, vcf, fai, seq_dict, indexed_ref_fai->
        tuple(fastq_id, bam, fasta, vcf, fai, seq_dict, indexed_ref_fai)
        }
    
    recal_file = createRecalModel(recalibration_inputs)

    recal_bam = recalBam(recal_file)

    vcf_file = callVariants(recal_bam)

    thousand_genome_vcf = downloadVariantArrays()

    indexed_vcf = indexCalledVariants(vcf_file)

    snp_vcf = purgeIndels(indexed_vcf)

    indexed_snp_vcf = indexSNPS(snp_vcf)

    filtered_variants = filterVariants(indexed_snp_vcf) 

    indexed_filtered_vcf = indexVCF(filtered_variants)

    refining_input = indexed_filtered_vcf
                    .combine(thousand_genome_vcf.variant_array_vcf.collect())
                    .combine(thousand_genome_vcf.variant_array_idx.collect())

    refined_vcf = calcGenotypePosteriors(refining_input)

    annotation_results = variantAnnotation(refined_vcf)
}



