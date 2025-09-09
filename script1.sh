//==============================================================================
// NEXTFLOW PIPELINE FOR BREAST CANCER SUSCEPTIBILITY GENE VARIANT ANALYSIS
// Adapted from GATK Best Practices for targeted exome sequencing analysis
//==============================================================================

//------------------------------------------------------------------------------
// PIPELINE PARAMETERS CONFIGURATION
//------------------------------------------------------------------------------
// Reference genome file (GRCh38 assembly, slim version for reduced size)
params.fasta_file = 'Homo_sapiens.GRCh38.dna_sm.slim.fa'

// Input sample identifier - currently configured for single sample analysis
params.fastq_id = ["ERR166333"]

// Base directory for reference files and output storage
params.index_dir = '/scratch/c.c24053373/test'

// Input VCF file path for variant filtering process
params.filtered_vcf_file = 'ERR166333_1.trimmed.fastq.sorted.markdup.rg.recal.vcf.gz'

// Population reference VCF for genotype posterior calculation (1000 Genomes data)
params.supporting_vcf = '1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf'

// Final filtered VCF file path for downstream processing
params.final_vcf = 'ERR166333_1.trimmed.fastq.sorted.markdup.rg.recal.snp.filtered.vcf.gz'

// Known variants VCF for base quality recalibration (Mills and 1000G gold standard)
params.vcf_file = 'Mills_and_1000G_gold_standard.indels.hg38.vcf'

// ANNOVAR database directory for variant functional annotation
params.annovar_db = '/scratch/c.c24053373/test/annovar/annovar/humandb'

//==============================================================================
// PROCESS DEFINITIONS - QUALITY CONTROL AND PREPROCESSING
//==============================================================================

/**
 * PROCESS: fastqTrim
 * PURPOSE: Remove low-quality bases and adapter sequences from raw FASTQ files
 * TOOL: Trimmomatic v0.39
 * INPUT: Paired-end FASTQ files (R1 and R2)
 * OUTPUT: Quality-trimmed FASTQ files (paired and unpaired)
 */
process fastqTrim{
    input:
    tuple path(fastq_1), path (fastq_2)  // Paired-end FASTQ input files

    output:
    path "*.trimmed.fastq.gz"            // Trimmed FASTQ output files

    script:
    """
    # Load required Java environment for Trimmomatic
    module load java/1.8

    # Trimmomatic paired-end trimming with quality and adapter filtering
    # PE: Paired-end mode
    # -phred33: Quality score encoding format
    # ILLUMINACLIP: Remove Illumina adapter sequences (TruSeq adapters)
    # LEADING:3: Remove leading low quality bases (below quality 3)
    # TRAILING:3: Remove trailing low quality bases (below quality 3)
    # SLIDINGWINDOW:4:15: Sliding window trimming (window size 4, avg quality 15)
    # MINLEN:36: Drop reads shorter than 36 bases after trimming
    java -jar /apps/genomics/trimmomatic/0.39/trimmomatic-0.39.jar PE -phred33 \
    ${fastq_1} ${fastq_2} \
    ${fastq_1.baseName}.trimmed.fastq.gz ${fastq_1.baseName}.unpaired.fastq.gz \
    ${fastq_2.baseName}.trimmed.fastq.gz ${fastq_2.baseName}.unpaired.fastq.gz \
    ILLUMINACLIP:/scratch/c.c24053373/test/TruSeq-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

/**
 * PROCESS: FASTQC
 * PURPOSE: Generate quality control reports for sequencing data
 * TOOL: FastQC v0.11.8
 * INPUT: Trimmed FASTQ files
 * OUTPUT: Quality control reports (HTML and ZIP formats)
 */
process FASTQC{
    // Publish results to output directory for review
    publishDir params.index_dir, mode: 'copy'

    input:
    path fastq_files                     // Trimmed FASTQ files for QC analysis

    output:
    path "*fastqc.zip"                   // Compressed QC data
    path "*fastqc.html"                  // Human-readable QC reports

    script:
    """
    # Load FastQC module
    module load FastQC/0.11.8
    
    # Generate quality control reports for all input FASTQ files
    fastqc ${fastq_files.join(' ')}
    """
}

//==============================================================================
// PROCESS DEFINITIONS - READ ALIGNMENT AND PROCESSING
//==============================================================================

/**
 * PROCESS: genomeMap
 * PURPOSE: Align trimmed reads to reference genome
 * TOOL: BWA-MEM v0.7.17
 * INPUT: Reference genome directory, FASTA filename, trimmed paired-end reads
 * OUTPUT: SAM alignment file
 */
process genomeMap{
    input:
    path data_dir                        // Directory containing reference genome
    val fasta_file                       // Reference genome filename
    tuple path(trimmed_fastq_1), path(trimmed_fastq_2)  // Trimmed paired-end reads

    output:
    path "*.sam"                         // SAM alignment output

    script:
    """
    # Load BWA alignment tool
    module load bwa/0.7.17
    
    # BWA-MEM alignment: Map paired-end reads to reference genome
    # mem: BWA-MEM algorithm (optimized for 70-100bp reads)
    # Output redirected to SAM format file
    bwa mem ${data_dir}/${fasta_file} ${trimmed_fastq_1} ${trimmed_fastq_2} > ${trimmed_fastq_1.baseName}.sam
    """
}

/**
 * PROCESS: sortBam
 * PURPOSE: Sort SAM file by genomic coordinates and convert to BAM format
 * TOOL: SAMtools v1.17
 * INPUT: SAM alignment file
 * OUTPUT: Coordinate-sorted BAM file
 */
process sortBam{
    input:
    path sam_file                        // Input SAM alignment file

    output:
    path "*.sorted.bam"                  // Coordinate-sorted BAM output

    script:
    """
    # Load SAMtools for BAM manipulation
    module load samtools/1.17
    
    # Sort SAM file by genomic coordinates and convert to compressed BAM format
    # -o: Output filename
    samtools sort -o ${sam_file.baseName}.sorted.bam ${sam_file}
    """
}

/**
 * PROCESS: markDuplicates
 * PURPOSE: Identify and mark PCR/optical duplicate reads
 * TOOL: Picard MarkDuplicates v2.27.5
 * INPUT: Coordinate-sorted BAM file
 * OUTPUT: BAM file with duplicates marked
 * NOTE: Duplicates are marked but not removed (REMOVE_DUPLICATES=false)
 */
process markDuplicates{
    input:
    path sorted_bam_file                 // Sorted BAM input

    output:
    path "*.sorted.markdup.bam"          // BAM with duplicates marked

    script:
    """
    # Load Java and Picard tools
    module load java/1.8
    module load picard/2.27.5
    
    # Mark duplicate reads using Picard MarkDuplicates
    # -Xmx32g: Allocate 32GB memory (essential for large datasets)
    # I: Input BAM file
    # O: Output BAM file with duplicates marked
    # M: Metrics file with duplicate statistics
    # REMOVE_DUPLICATES=false: Mark but don't remove duplicates
    java -Xmx32g -jar \$PICARD MarkDuplicates \
        I=${sorted_bam_file} \
        O=${sorted_bam_file.baseName}.markdup.bam \
        M=${sorted_bam_file.baseName}.metrics \
        REMOVE_DUPLICATES=false
    """
}

/**
 * PROCESS: bamStats
 * PURPOSE: Generate alignment statistics for quality assessment
 * TOOL: BAMtools v170119
 * INPUT: BAM file with marked duplicates
 * OUTPUT: Alignment statistics file
 */
process bamStats {
    input:
    path markdup_bam_file                // BAM file with marked duplicates

    output:
    path "*.sorted.markdup.stats"        // Alignment statistics output

    script:
    """
    # Load BAMtools for alignment statistics
    module load bamtools/170119
    
    # Generate comprehensive alignment statistics
    # Including: read counts, mapping rates, insert sizes, etc.
    bamtools stats -in ${markdup_bam_file} > ${markdup_bam_file.baseName}.stats
    """
}

/**
 * PROCESS: indexBam
 * PURPOSE: Create BAM index for random access (required for GATK tools)
 * TOOL: SAMtools v1.17
 * INPUT: BAM file with marked duplicates
 * OUTPUT: BAM index file (.bai)
 */
process indexBam {
    input:
    path markdup_bam_file                // BAM file to index

    output:
    path "*.sorted.markdup.bam.bai"      // BAM index output

    script:
    """
    # Load SAMtools for indexing
    module load samtools/1.17
    
    # Create BAM index for random access
    # Required by GATK tools for efficient BAM file access
    samtools index ${markdup_bam_file}
    """
}

/**
 * PROCESS: addreadGroups
 * PURPOSE: Add read group metadata required by GATK tools
 * TOOL: Picard AddOrReplaceReadGroups v2.27.5
 * INPUT: BAM file with marked duplicates
 * OUTPUT: BAM file with read group information
 */
process addreadGroups {
    input:
    path markdup_bam_file                // Input BAM file

    output:
    path "*.sorted.markdup.rg.bam"       // BAM with read groups added

    script:
    """
    # Load Picard and SAMtools
    module load picard/2.27.5
    module load samtools/1.17
    
    # Add read group metadata required by GATK
    # SO=coordinate: Sort order is coordinate-based
    # RGID: Read group identifier
    # RGLB: Library identifier
    # RGPL: Platform (illumina)
    # RGSM: Sample name (ERR166333)
    # RGPU: Platform unit identifier
    # CREATE_INDEX=true: Create BAM index automatically
    java -jar \$PICARD AddOrReplaceReadGroups \
        INPUT=${markdup_bam_file} \
        OUTPUT=${markdup_bam_file.baseName}.rg.bam \
        SO=coordinate \
        RGID=1 \
        RGLB=lib1 \
        RGPL=illumina \
        RGSM=ERR166333 \
        RGPU=unit1 \
        CREATE_INDEX=true
    """
}

//==============================================================================
// PROCESS DEFINITIONS - REFERENCE GENOME PREPARATION
//==============================================================================

/**
 * PROCESS: indexReference
 * PURPOSE: Create FASTA index for reference genome (required for GATK)
 * TOOL: SAMtools v1.3.1
 * INPUT: Reference FASTA file
 * OUTPUT: FASTA index (.fai file)
 */
process indexReference {
    // Publish index to output directory for reuse
    publishDir params.index_dir, mode: 'copy'

    input: 
    path fasta_file                      // Reference FASTA file

    output:
    path "*.fai"                         // FASTA index output

    script:
    """
    # Load SAMtools for FASTA indexing
    module load samtools/1.3.1
    
    # Create FASTA index for random access to reference genome
    # Required by GATK tools for efficient reference access
    samtools faidx ${fasta_file}
    """
}

/**
 * PROCESS: indexKnownVariants
 * PURPOSE: Create index for known variants VCF file (required for BQSR)
 * TOOL: GATK v4.1.2.0
 * INPUT: Known variants VCF file
 * OUTPUT: VCF index (.idx file)
 */
process indexKnownVariants {
    // Publish index to output directory for reuse
    publishDir params.index_dir, mode: 'copy'

    input:
    path vcf_file                        // Known variants VCF file

    output:
    path "${vcf_file}.idx"               // VCF index output

    script:
    """
    # Load GATK for VCF indexing
    module load GATK/4.1.2.0
    
    # Create VCF index for known variants file
    # Required for Base Quality Score Recalibration (BQSR)
    gatk IndexFeatureFile -F ${vcf_file}
    """
}

/**
 * PROCESS: createSequenceDictionary
 * PURPOSE: Create sequence dictionary for reference genome (required by GATK)
 * TOOL: Picard CreateSequenceDictionary v2.27.5
 * INPUT: Reference FASTA file
 * OUTPUT: Sequence dictionary (.dict file)
 */
process createSequenceDictionary {
    // Publish dictionary to output directory for reuse
    publishDir params.index_dir, mode: 'copy'

    input: 
    path fasta_file                      // Reference FASTA file

    output: 
    path "*.dict"                        // Sequence dictionary output

    script:
    """
    # Load Picard for sequence dictionary creation
    module load picard/2.27.5
    
    # Create sequence dictionary containing chromosome information
    # R: Reference FASTA file
    # O: Output dictionary file
    # Required by GATK tools for reference genome validation
    java -jar \$PICARD CreateSequenceDictionary \
        R=${fasta_file} \
        O=${fasta_file.baseName}.dict
    """
}

//==============================================================================
// PROCESS DEFINITIONS - BASE QUALITY SCORE RECALIBRATION (BQSR)
//==============================================================================

/**
 * PROCESS: baseRecalibrator
 * PURPOSE: Analyze systematic base quality score errors for recalibration
 * TOOL: GATK BaseRecalibrator v4.1.2.0
 * INPUT: BAM with read groups, reference genome, known variants
 * OUTPUT: Recalibration table with quality score corrections
 */
process baseRecalibrator {
    input: 
    path bam_file                        // BAM file with read groups
    path fasta_file                      // Reference genome filename
    path known_sites_vcf                 // Known variants for recalibration

    output:
    path "*.recal_data.table"            // Base recalibration table

    script:
    """
    # Load GATK for base quality recalibration
    module load GATK/4.1.2.0
    
    # Analyze systematic base quality score errors
    # -I: Input BAM file
    # -R: Reference genome
    # --known-sites: Known variants to exclude from analysis
    # -O: Output recalibration table
    gatk BaseRecalibrator \
        -I ${bam_file} \
        -R ${params.index_dir}/${fasta_file} \
        --known-sites ${params.index_dir}/${known_sites_vcf} \
        -O ${bam_file.baseName}.recal_data.table
    """
}

/**
 * PROCESS: applyBQSR
 * PURPOSE: Apply base quality score recalibration to BAM file
 * TOOL: GATK ApplyBQSR v4.1.2.0
 * INPUT: BAM file, reference genome, recalibration table
 * OUTPUT: BAM file with recalibrated base quality scores
 */
process applyBQSR{
    input:
    path bam_file                        // Input BAM file
    path fasta_file                      // Reference genome filename
    path recal_table                     // Base recalibration table

    output:
    path "*.recal.bam"                   // BAM with recalibrated qualities

    script:
    """
    # Load GATK for applying base quality recalibration
    module load GATK/4.1.2.0
    
    # Apply base quality score recalibration
    # -R: Reference genome
    # -I: Input BAM file
    # --bqsr-recal-file: Recalibration table from BaseRecalibrator
    # -O: Output BAM with corrected base qualities
    gatk ApplyBQSR \
        -R ${params.index_dir}/${fasta_file} \
        -I ${bam_file} \
        --bqsr-recal-file ${recal_table} \
        -O ${bam_file.baseName}.recal.bam
    """
}

//==============================================================================
// PROCESS DEFINITIONS - VARIANT CALLING AND INDEXING
//==============================================================================

/**
 * PROCESS: haplotypeCaller
 * PURPOSE: Call SNPs and INDELs using GATK's local assembly approach
 * TOOL: GATK HaplotypeCaller v4.1.2.0
 * INPUT: Recalibrated BAM file, reference genome
 * OUTPUT: VCF file with variant calls
 */
process haplotypeCaller {
    // Publish VCF to output directory
    publishDir params.index_dir, mode: 'copy'

    input:
    path bam_file                        // Recalibrated BAM file
    path fasta_file                      // Reference genome filename

    output:
    path "*.vcf.gz"                      // Compressed VCF output

    script:
    """
    # Load GATK and Java for variant calling
    module load GATK/4.1.2.0
    module load java/1.8
    
    # Call variants using HaplotypeCaller
    # --java-options "-Xmx4g": Allocate 4GB memory
    # -R: Reference genome
    # -I: Input recalibrated BAM file
    # -O: Output VCF file (compressed)
    gatk --java-options "-Xmx4g" HaplotypeCaller \
        -R ${params.index_dir}/${fasta_file} \
        -I ${bam_file} \
        -O ${bam_file.baseName}.vcf.gz
    """
}

/**
 * PROCESS: indexVCF
 * PURPOSE: Create index for VCF file (required for downstream tools)
 * TOOL: GATK IndexFeatureFile v4.1.2.0
 * INPUT: VCF file from HaplotypeCaller
 * OUTPUT: VCF index (.tbi file)
 */
process indexVCF {
    // Publish index to output directory
    publishDir params.index_dir, mode: 'copy'

    input: 
    path vcf_file                        // VCF file to index

    output:
    path "*.vcf.gz.tbi"                  // VCF index output

    script:
    """
    # Load GATK and Java for VCF indexing
    module load GATK/4.1.2.0
    module load java/1.8

    # Create VCF index for random access
    # Required by downstream GATK tools
    gatk IndexFeatureFile -F ${vcf_file}
    """
}

//==============================================================================
// PROCESS DEFINITIONS - VARIANT FILTERING AND REFINEMENT
//==============================================================================

/**
 * PROCESS: filterVariants
 * PURPOSE: Select SNPs and apply hard quality filters
 * TOOL: GATK SelectVariants + VariantFiltration v4.1.2.0
 * INPUT: Trigger value (workflow coordination)
 * OUTPUT: Filtered SNP VCF file
 */
process filterVariants {
    // Publish filtered VCF to output directory
    publishDir params.index_dir, mode: 'copy'

    input:
    val trigger                          // Dummy input for workflow coordination

    output:
    path "*.snp.filtered.vcf.gz"         // Filtered SNP VCF output

    script:
    // Extract base filename for consistent naming
    def filename = params.filtered_vcf_file
    def base_name = filename.substring(0, filename.lastIndexOf('.vcf.gz'))

    """
    # Load GATK and Java for variant filtering
    module load GATK/4.1.2.0
    module load java/1.8

    # Step 1: Select only SNPs from the VCF file
    # --variant: Input VCF file with all variants
    # --select-type SNP: Select only single nucleotide polymorphisms
    # --output: Output VCF with SNPs only
    gatk SelectVariants \
        --variant ${params.index_dir}/${params.filtered_vcf_file} \
        --select-type SNP \
        --output ${base_name}.snp.vcf.gz

    # Step 2: Apply hard quality filters to SNPs
    # Filter expressions based on GATK best practices for exome data:
    # QD<2.0: Quality by Depth (variant confidence/read depth)
    # QUAL<30.0: Overall variant quality score
    # SOR>3.0: Strand Odds Ratio (strand bias detection)
    # FS>60.0: Fisher Strand test (strand bias detection)
    # MQ<40.0: Mapping Quality (alignment confidence)
    # MQRankSum<-12.5: Mapping Quality Rank Sum test
    # ReadPosRankSum<-8.0: Read Position Rank Sum test
    gatk VariantFiltration \
        --variant ${base_name}.snp.vcf.gz \
        --filter-expression "QD<2.0" --filter-name "QD2" \
        --filter-expression "QUAL<30.0" --filter-name "QUAL30" \
        --filter-expression "SOR>3.0" --filter-name "SOR3" \
        --filter-expression "FS>60.0" --filter-name "FS60" \
        --filter-expression "MQ<40.0" --filter-name "MQ40" \
        --filter-expression "MQRankSum<-12.5" --filter-name "MQRankSum-12.5" \
        --filter-expression "ReadPosRankSum<-8.0" --filter-name "ReadPosRankSum-8" \
        --output ${base_name}.snp.filtered.vcf.gz
    """
}

/**
 * PROCESS: indexFilteredVCF
 * PURPOSE: Create index for filtered VCF file
 * TOOL: GATK IndexFeatureFile v4.1.2.0
 * INPUT: Trigger value (workflow coordination)
 * OUTPUT: Index for filtered VCF file
 */
process indexFilteredVCF {
    // Publish index to output directory
    publishDir params.index_dir, mode: 'copy'

    input:
    val trigger                          // Dummy input for workflow coordination

    output:
    path "*.vcf.gz.tbi"                  // Filtered VCF index output

    script:
    // Extract base filename for consistent naming
    def filename = params.final_vcf
    def base_name = filename.substring(0, filename.lastIndexOf('.vcf.gz'))

    """
    # Load GATK and Java for VCF indexing
    module load GATK/4.1.2.0
    module load java/1.8

    # Create symbolic link to filtered VCF in working directory
    # Required because GATK creates index in same directory as VCF
    ln -sf ${params.index_dir}/${params.final_vcf} ${base_name}.vcf.gz

    # Create index for filtered VCF file
    gatk IndexFeatureFile -F ${base_name}.vcf.gz
    """
}

/**
 * PROCESS: calculateGenotypePosteriors
 * PURPOSE: Refine genotype calls using population allele frequencies
 * TOOL: GATK CalculateGenotypePosteriors v4.1.2.0
 * INPUT: Trigger value, population reference VCF and index
 * OUTPUT: VCF with refined genotype posterior probabilities
 */
process calculateGenotypePosteriors {
    // Publish refined VCF to output directory
    publishDir params.index_dir, mode: 'copy'
    
    input:
    val trigger                          // Dummy input for workflow coordination
    path supporting_vcf                  // Population reference VCF (1000 Genomes)
    path supporting_vcf_idx              // Population reference VCF index

    output:
    path "*.ref.vcf.gz"                  // VCF with posterior probabilities

    script:
    // Extract base filename for consistent naming
    def filename = params.final_vcf
    def base_name = filename.substring(0, filename.lastIndexOf('.vcf.gz'))

    """
    # Load GATK and Java for posterior calculation
    module load GATK/4.1.2.0
    module load java/1.8

    # Calculate genotype posterior probabilities using population data
    # --java-options "-Xmx4g": Allocate 4GB memory
    # -V: Input filtered VCF file
    # -O: Output VCF with posterior probabilities
    # -supporting: Population reference VCF (1000 Genomes Project)
    gatk --java-options "-Xmx4g" CalculateGenotypePosteriors \
        -V ${params.index_dir}/${params.final_vcf} \
        -O ${base_name}.ref.vcf.gz \
        -supporting ${supporting_vcf}
    """
}

//==============================================================================
// PROCESS DEFINITIONS - VARIANT ANNOTATION
//==============================================================================

/**
 * PROCESS: annotateVariants
 * PURPOSE: Add functional annotations to variants using multiple databases
 * TOOL: ANNOVAR v20210525
 * INPUT: Trigger value, ANNOVAR database directory
 * OUTPUT: Comprehensively annotated VCF file
 */
process annotateVariants {
    // Publish annotated VCF to output directory
    publishDir params.index_dir, mode: 'copy'

    input:
    val trigger                          // Dummy input for workflow coordination
    path annovar_db                      // ANNOVAR database directory

    output:
    path "*.annot.vcf.gz"                // Annotated VCF output

    script:
    // Extract base filename for consistent naming
    def filename = params.final_vcf
    def base_name = filename.substring(0, filename.lastIndexOf('.vcf.gz'))

    """
    # Load ANNOVAR for variant annotation
    module load annovar/20210525

    # Annotate variants using multiple ANNOVAR databases
    # Input: VCF with refined genotype posteriors
    # -buildver hg38: Use GRCh38/hg38 genome build
    # -out: Output file prefix
    # -remove: Remove temporary files after completion
    # -protocol: Annotation databases to use:
    #   - refGene: Gene-based annotation (transcripts, consequences)
    #   - cytoBand: Cytogenetic band information
    #   - avsnp150: dbSNP build 150 (variant IDs, frequencies)
    #   - clinvar_20210501: ClinVar clinical significance (May 2021)
    # -operation: Annotation types (g=gene-based, r=region-based, f=filter-based)
    # -nastring .: Use '.' for missing annotations
    # -vcfinput: Input is in VCF format
    table_annovar.pl \
        ${params.index_dir}/${base_name}.ref.vcf.gz \
        ${annovar_db} \
        -buildver hg38 \
        -out ${base_name}.ref \
        -remove \
        -protocol refGene,cytoBand,avsnp150,clinvar_20210501 \
        -operation g,r,f,f \
        -nastring . \
        -vcfinput

    # Compress the annotated VCF file and rename for consistency
    gzip ${base_name}.ref.hg38_multianno.vcf
    mv ${base_name}.ref.hg38_multianno.vcf.gz ${base_name}.ref.annot.vcf.gz
    """
}

//==============================================================================
// WORKFLOW DEFINITION - PROCESS ORCHESTRATION AND DATA FLOW
//==============================================================================

workflow {
    //--------------------------------------------------------------------------
    // INPUT DATA PREPARATION
    //--------------------------------------------------------------------------
    // Create channel for FASTQ input files
    // Maps sample IDs to paired-end FASTQ file paths
    fastq_ch = Channel
               .from(params.fastq_id)
               // Map sample ID to paired FASTQ file paths
               .map { id -> tuple("/scratch/c.c24053373/test/${id}_1.fastq", "/scratch/c.c24053373/test/${id}_2.fastq") }

    //--------------------------------------------------------------------------
    // PHASE 1: QUALITY CONTROL AND READ PREPROCESSING
    //--------------------------------------------------------------------------
    // Trim low-quality bases and adapters from raw FASTQ files
    trimmed_fastq_ch = fastqTrim(fastq_ch)

    // Generate quality control reports for trimmed reads
    FASTQC(trimmed_fastq_ch)

    //--------------------------------------------------------------------------
    // PHASE 2: READ ALIGNMENT AND BAM PROCESSING
    //--------------------------------------------------------------------------
    // Set up channels for reference genome alignment
    index_ch = Channel.fromPath(params.index_dir)        // Reference directory
    ref_ch = Channel.of(params.fasta_file)               // Reference filename
    
    // Align trimmed reads to reference genome using BWA-MEM
    sam_ch = genomeMap(index_ch, ref_ch, trimmed_fastq_ch)

    // Convert SAM to sorted BAM format
    sorted_bam_ch = sortBam(sam_ch)

    // Mark PCR and optical duplicate reads
    markdup_bam_ch = markDuplicates(sorted_bam_ch)

    // Generate alignment statistics for quality assessment
    bamStats(markdup_bam_ch)

    // Create BAM index for random access (required for IGV visualization)
    bamIndex_ch = indexBam(markdup_bam_ch)

    // Add read group metadata required by GATK tools
    readgroup_bam_ch = addreadGroups(markdup_bam_ch)

    //--------------------------------------------------------------------------
    // PHASE 3: REFERENCE GENOME PREPARATION (PARALLEL PROCESSING)
    //--------------------------------------------------------------------------
    // These processes run in parallel as they're independent
    
    // Create FASTA index for reference genome
    fasta_index_ch = indexReference(Channel.fromPath("${params.index_dir}/${params.fasta_file}"))

    // Create index for known variants VCF (required for BQSR)
    known_variants_ch = indexKnownVariants(Channel.fromPath("${params.index_dir}/${params.vcf_file}"))

    // Create sequence dictionary for reference genome
    fasta_dict_ch = createSequenceDictionary(Channel.fromPath("${params.index_dir}/${params.fasta_file}"))

    //--------------------------------------------------------------------------
    // PHASE 4: BASE QUALITY SCORE RECALIBRATION (BQSR)
    //--------------------------------------------------------------------------
    // Analyze systematic base quality score errors
    recal_table_ch = baseRecalibrator(
        readgroup_bam_ch,                                // BAM with read groups
        Channel.fromPath("${params.fasta_file}"),        // Reference genome
        Channel.fromPath("${params.vcf_file}")           // Known variants
    )

    // Apply base quality score recalibration to BAM file
    recal_bam_ch = applyBQSR(
        readgroup_bam_ch,                                // Input BAM
        Channel.fromPath("${params.fasta_file}"),        // Reference genome
        recal_table_ch                                   // Recalibration table
    )

    //--------------------------------------------------------------------------
    // PHASE 5: VARIANT CALLING AND INDEXING
    //--------------------------------------------------------------------------
    // Call variants using GATK HaplotypeCaller
    vcf_ch = haplotypeCaller(
        recal_bam_ch,                                    // Recalibrated BAM
        Channel.fromPath("${params.fasta_file}")         // Reference genome
    )
    
    // Create index for raw VCF file
    indexed_vcf_ch = indexVCF(vcf_ch)

    //--------------------------------------------------------------------------
    // PHASE 6: VARIANT FILTERING AND REFINEMENT (SEQUENTIAL PROCESSING)
    //--------------------------------------------------------------------------
    // Create workflow coordination triggers to ensure proper execution order
    
    // Wait for variant calling to complete before filtering
    haplotype_done = vcf_ch.collect().map { "step1_done" }
    
    // Filter variants (select SNPs and apply quality filters)
    filtered_vcf_ch = filterVariants(haplotype_done)
    
    // Wait for variant filtering to complete before indexing
    filter_done = filtered_vcf_ch.collect().map { "step2_done" }
    
    // Create index for filtered VCF file
    indexed_filtered_vcf_ch = indexFilteredVCF(filter_done)

    // Wait for filtered VCF indexing to complete
    index_done = indexed_filtered_vcf_ch.collect().map { "step3_done" }

    //--------------------------------------------------------------------------
    // PHASE 7: GENOTYPE REFINEMENT USING POPULATION DATA
    //--------------------------------------------------------------------------
    // Calculate genotype posterior probabilities using 1000 Genomes data
    posterior_vcf_ch = calculateGenotypePosteriors(
        index_done,                                      // Wait for previous step
        Channel.fromPath("${params.index_dir}/${params.supporting_vcf}"),     // 1000G VCF
        Channel.fromPath("${params.index_dir}/${params.supporting_vcf}.idx")  // 1000G index
    )
    
    //--------------------------------------------------------------------------
    // PHASE 8: COMPREHENSIVE VARIANT ANNOTATION
    //--------------------------------------------------------------------------
    // Wait for genotype posterior calculation to complete
    posterior_done = posterior_vcf_ch.collect().map { "step4_done" }
    
    // Add functional annotations using ANNOVAR databases
    annotated_vcf_ch = annotateVariants(
        posterior_done,                                  // Wait for previous step
        Channel.fromPath("${params.annovar_db}")         // ANNOVAR database directory
    )
}

//==============================================================================
// END OF PIPELINE
// Final output: Comprehensively annotated VCF file ready for clinical interpretation
// and breast cancer susceptibility gene analysis
//==============================================================================
