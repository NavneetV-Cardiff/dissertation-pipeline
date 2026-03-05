#!/bin/bash
#SBATCH --job-name=nextflow_pipeline
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=8:00:00
#SBATCH --output=nextflow_%j.log
#SBATCH --mail-user=(your_email_here)
#SBATCH --mail-type=ALL

module load Nextflow/25.10.2
nextflow run exomeseq.nf 
