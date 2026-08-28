#!/bin/bash
#TR14_QC_script.sh
#NOTE: this is just one sample as an example

module load fastqc
module load trimmomatic
cd /workspace/rweed/TR_phage_redo
#pre-QC assessment:
fastqc 00_SAMPLES/TR14-B05-S4_1.fastq 00_SAMPLES/TR14-B05-S4_2.fastq -o 01_QC/pre_QC_reports/
fastqc 00_SAMPLES/TR14-B05-S5_1.fastq 00_SAMPLES/TR14-B05-S5_2.fastq -o 01_QC/pre_QC_reports/

#trimmomatic
trimmomatic PE -trimlog TR14_B05_S4_trim.log -phred33 00_SAMPLES/TR14-B05-S4_1.fastq 00_SAMPLES/TR14-B05-S4_2.fastq -baseout 01_QC/TR14_B05_S4-QUALITY_PASSED_trimmed ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75
trimmomatic PE -trimlog TR14_B05_S5_trim.log -phred33 00_SAMPLES/TR14-B05-S5_1.fastq 00_SAMPLES/TR14-B05-S5_2.fastq -baseout 01_QC/TR14_B05_S5-QUALITY_PASSED_trimmed ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75

#post-QC assessment:
fastqc 01_QC/TR14_B05_S4-QUALITY_PASSED_trimmed_1P 01_QC/TR14_B05_S4-QUALITY_PASSED_trimmed_2P -o 01_QC/post_QC_reports/
fastqc 01_QC/TR14_B05_S5-QUALITY_PASSED_trimmed_1P 01_QC/TR14_B05_S5-QUALITY_PASSED_trimmed_2P -o 01_QC/post_QC_reports/
