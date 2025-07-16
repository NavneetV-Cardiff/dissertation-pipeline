# Load SRAtoolkit
module load sra-tools/2.10.0 
  # Download SRA data
  prefetch ERR166333
  fasterq-dump ERR166333
  # prefetch ERR166337
  # fasterq-dump ERR166337

  # scp ERR166333 -S 

 gzip ERR166333*fastq 


# Download the Human Genome
 wget "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"  
  gunzip Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz 

# Download the cancer genomes

 wget "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa.gz" 

 wget "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.3.fa.gz" 

 wget "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.4.fa.gz" 

 wget "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.5.fa.gz" 

 wget "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.9.fa.gz" 

 wget "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.14.fa.gz" 

 wget "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.16.fa.gz" 

 wget "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.17.fa.gz" 

 wget "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.22.fa.gz" 

 gunzip *.gz 

 cat *chromosome*fa > Homo_sapiens.GRCh38.dna_sm.slim.fa 


# Indexing the FASTA reference genome before mapping with BWA

# Just checking the .slim file to see which chromosomes have been included
grep '>' Homo_sapiens.GRCh38.dna_sm.slim.fa 

#### STEP 1

module load bwa/0.7.17 

# This should def be part of the final pipeline
# Indexing FASTA file
bwa index -a bwtsw Homo_sapiens.GRCh38.dna_sm.slim.fa 

####

#### STEP 2
# Trim the FASTQ files
module load java/1.8 
module load trimmomatic/0.39 

# Here the fastq files are test files which are "slim"
java -jar /apps/genomics/trimmomatic/0.39/trimmomatic-0.39.jar PE -phred33 ERR166333.slim_1.fastq.gz ERR166333.slim_2.fastq.gz ERR166333.slim.trim_1.fastq.gz ERR166333.slim.trim.unpaired_1.fastq.gz ERR166333.slim.trim_2.fastq.gz ERR166333.slim.trim.unpaired_2.fastq.gz ILLUMINACLIP:TruSeq-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 

####

# Load nextflow
module load nextflow/24.10.4

#### STEP 3
# QC of FASTQ files
module load FastQC/0.11.8 

fastqc ERR166333.slim_1.fastq.gz ERR166333.slim_2.fastq.gz ERR166333.slim.trim_1.fastq.gz ERR166333.slim.trim_2.fastq.gz 
#### 

#### STEP 4
# Mapping FASTQ reads to reference sequence

module load bwa/0.7.17 
bwa mem Homo_sapiens.GRCh38.dna_sm.slim.fa ERR166333.slim.trim_1.fastq.gz ERR166333.slim.trim_2.fastq.gz > ERR166333.slim.sam 

####


#### STEP 5
# Creating sorted BAM file

module load samtools/1.17 
samtools view -bS ERR166333.slim.sam > ERR166333.slim.bam 
# Check difference in sizes
ls -l ERR166333.slim.*am 

samtools sort -o ERR166333.slim.sorted.bam ERR166333.slim.bam 

####

#### STEP 6
# Marking Duplicates
module load java/1.8 

module load picard/2.27.5

java -jar $PICARD MarkDuplicates I=ERR166333.slim.sorted.bam O=ERR166333.slim.sorted.markdup.bam M=ERR166333.slim.sorted.metrics REMOVE_DUPLICATES=false

####

### STEP 7
# Second QC (on BAMtools)

module load bamtools/170119 
bamtools stats -in ERR166333.slim.sorted.markdup.bam > ERR166333.slim.sorted.markdup.stats

###

#### STEP 8
# Creating an Indexed BAM and opening in IGV
module load samtools/1.17 
samtools index ERR166333.slim.sorted.markdup.bam


#####

#### STEP 9
# adding read groups to the BAM file
module load picard/2.27.5 
module load samtools/1.17 
java -jar $PICARD AddOrReplaceReadGroups \
INPUT=ERR166333.slim.sorted.markdup.bam \
OUTPUT=ERR166333.slim.sorted.markdup.rg.bam SO=coordinate RGID=1 RGLB=lib1 RGPL=illumina RGSM=ERR166333 RGPU=unit1 CREATE_INDEX=True 

# Not needed in final script
cp ERR166333.slim.sorted.markdup.rg.bam temp.rg.bam 
samtools view temp.rg.bam > temp.rg.sam 
more temp.rg.sam 


#####

#### STEP 10
# PHRED score recalibration
module load samtools/1.17 

module load GATK/4.1.2.0 

module load picard/2.27.5

samtools faidx Homo_sapiens.GRCh38.dna_sm.slim.fa 

 

ls Homo_sapiens.GRCh38.dna_sm.slim.fa.fai

wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz 

# Index known variants
gatk IndexFeatureFile -F Mills_and_1000G_gold_standard.indels.hg38.vcf.gz 

# Create a sequence dictionary from the fasta file
java -jar $PICARD CreateSequenceDictionary R=Homo_sapiens.GRCh38.dna_sm.slim.fa O=Homo_sapiens.GRCh38.dna_sm.slim.dict 

more Homo_sapiens.GRCh38.dna_sm.slim.dict

# Create the model for the recalibration step 
gatk BaseRecalibrator -I ERR166333.slim.sorted.markdup.rg.bam -R Homo_sapiens.GRCh38.dna_sm.slim.fa --known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -O ERR166333.recal_data.table 

# Apply the model to the newly recalibrated file
gatk ApplyBQSR -R Homo_sapiens.GRCh38.dna_sm.slim.fa -I ERR166333.slim.sorted.markdup.rg.bam --bqsr-recal-file ERR166333.recal_data.table -O ERR166333.slim.sorted.markdup.rg.recal.bam



####


#### STEP 11
# Calling variants
module load GATK/4.1.2.0 

module load java/1.8

#
gatk --java-options "-Xmx4g" HaplotypeCaller -R Homo_sapiens.GRCh38.dna_sm.slim.fa -I ERR166333.slim.sorted.markdup.rg.recal.bam -O ERR166333.slim.vcf.gz 

ls ERR166333.slim.vcf.gz 



####


#### Step 12
# Filter variants
module load GATK/4.1.2.0 

module load java/1.8 
gatk SelectVariants --variant ERR166333.slim.vcf.gz --select-type SNP --output ERR166333.slim.snp.vcf.gz 
gatk VariantFiltration  --variant ERR166333.slim.snp.vcf.gz --filter-expression "QD<2.0" --filter-name "QD2" --filter-expression "QUAL<30.0" --filter-name "QUAL30"  --filter-expression "SOR>3.0" --filter-name "SOR3" --filter-expression "FS>60.0" --filter-name "FS60" --filter-expression "MQ<40.0" --filter-name "MQ40" --filter-expression "MQRankSum<-12.5" --filter-name "MQRankSum-12.5" --filter-expression "ReadPosRankSum<-8.0" --filter-name "ReadPosRankSum-8" --output ERR166333.slim.snp.filtered.vcf.gz



####

#### STEP 13
#
module load GATK/4.1.2.0 

module load java/1.8 

# Download reference .idx files from 1000genomes project

wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf 

 

wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf.idx

# Calculate posterior probabilities
gatk --java-options "-Xmx4g" CalculateGenotypePosteriors -V ERR166333.slim.snp.filtered.vcf.gz -O ERR166333.slim.snp.filtered.ref.vcf.gz -supporting 1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf

####


#### STEP 14
# Annotate variants

gatk FuncotatorDataSourceDownloader --germline --validate-integrity --extract-after-download 

# Run funcotator to add gene annotations to vcf
gatk Funcotator --variant ERR166333.slim.snp.filtered.ref.vcf.gz --reference Homo_sapiens.GRCh38.dna_sm.slim.fa --ref-version hg38 --data-sources-path funcotator_dataSources.v1.6.20190124g --output ERR166333.slim.snp.filtered.ref.annot.vcf.gz --output-file-format VCF 



####





