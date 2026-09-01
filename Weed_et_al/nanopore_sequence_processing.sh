#!/bin/bash

##CONVERT TO POD5##
#path to fast5 files of run1
cd /Users/mollymoynihan/Documents/Nanopore/TR14_Oct2022/RUFF_TR_32_10272022_run1/no_sample/20221027_1554_MN40625_FAS63891_617524a9

# Note the different --one-to-one path which is now the current working directory.
# The new sub-directory output_pod5/input is created.
pod5 convert fast5 fast5/*.fast5 --output output_pod5s --one-to-one ./
ls output_pod5s/

#----------------------------------
#path to fast5 files of run2
cd /Users/mollymoynihan/Documents/Nanopore/TR14_Oct2022/RUFF_TR_32_10272022_run2/RUFF_TR_32_10272022/no_sample/20221028_1303_MN40625_FAS63891_f375616a

# Note the different --one-to-one path which is now the current working directory.
# The new sub-directory output_pod5/input is created.
pod5 convert fast5 fast5/*.fast5 --output output_pod5s --one-to-one ./
ls output_pod5s/

##BASECALLING##
cd /Users/mollymoynihan/Documents/Nanopore/TR14_Oct2022/pod5_files
# run1
/Users/mollymoynihan/Downloads/dorado-0.5.3-osx-arm64/bin/dorado basecaller hac run1/ > run1_calls.bam

# run2
/Users/mollymoynihan/Downloads/dorado-0.5.3-osx-arm64/bin/dorado basecaller hac run2/ > run2_calls.bam

# Convert bam files to fastq files 
samtools fastq run1_calls.bam > run1.fastq
samtools fastq run2_calls.bam > run2.fastq

# Concatenate files (file size is 24GB)
cat run1.fastq run2.fastq > TR14_all_runs.fastq

# Moved fastq files to a new folder called 'fastq_files' 

##QUALITY CONTROL##

cd /Users/mollymoynihan/Documents/Nanopore/TR14_Oct2022/fastq_files

filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 TR14_all_runs.fastq | gzip > TR14_all_QCfiltered.fastq.gz

##ASSEMBLY##
# Running Flye with standard settings, identifying the fastq files as "corrected" because I ran Fitlong above 
python3 Flye/bin/flye --nano-corr Nanopore/TR14_Oct2022/fastq_files/TR14_all_QCfiltered.fastq.gz --out-dir Nanopore/TR14_Oct2022/Flye_assembly --threads 6

##SEQUENCE CORRECTION##
module load medaka

cd /workspace/mmoynihan/Nanopore/TR14_Assembly

source ${MEDAKA}  # i.e. medaka/venv/bin/activate
NPROC=$(nproc)
BASECALLS=02_QC/TR14_all_QCfiltered.fastq.gz
DRAFT=03_ASSEMBLY/Flye_assembly/assembly.fasta
OUTDIR=04_LR_POLISHING/medaka_consensus
medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR} -t 20 -m r941_min_high_g303

##POLISHING##
 
#Copied consensus.fasta from LR polishing into new SR polishing folder 05_

cd 05_SR_POLISHING

# align the reads to your draft genome with BWA, use -a flag to align all reads to all possible locations
bwa index consensus.fasta
bwa mem -t 8 -a consensus.fasta Illumina/TR14-B05_1.fastq.gz > alignments_1.sam
bwa mem -t 8 -a consensus.fasta Illumina/TR14-B05_2.fastq.gz > alignments_2.sam

# Filter the alignments using Polypolish's insert size filter. This step is optional but recommended:
polypolish filter --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam

# Then give the draft genome and the alignments to Polypolish. It will output information to stderr and the polished assembly to stdout, so redirect its output to a file:
polypolish polish consensus.fasta filtered_1.sam filtered_2.sam > polished.fasta

# Clean up the index files (made by bwa index) and alignments to save disk space:
rm *.amb *.ann *.bwt *.pac *.sa *.sam
