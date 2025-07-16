# Downloading the data

# Load SRAtoolkit
module load sra-tools/2.10.0 
  # Download SRA data
  prefetch ERR166333
  fasterq-dump ERR166333

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

