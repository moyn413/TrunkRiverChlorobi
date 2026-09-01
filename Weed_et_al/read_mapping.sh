#!/bin/bash
##PREP MAGS##
#These MAGs have been quality filtered as outlined in host_matching.sh
#They have also been filtered for MAGs derived from B05 only
#dereplicate using dRep
cd /workspace/rweed/TR_phage_redo/06_HOSTS
/workspace/rweed/TR_phage_redo/UTILS/drep/bin/dRep dereplicate B05_bins_medium_and_up_dereplicated -g B05_bins_medium_and_up/*.fa
mkdir 10_MAPPING/MAGs_only_derep
cp 06_HOSTS/B05_bins_medium_and_up_dereplicated/dereplicated_genomes/*.fa 10_MAPPING/MAGs_only_derep/
#add contig_20
cp 11_NANOPORE/01_ASSEMBLY/polished_contigs_separated/contig_20.fasta 10_MAPPING/MAGs_only_derep/
#change all extensions to be .fasta
rename fa fasta *.fa

##PREP vMAGS##
#These vMAGs are the quality controlled, binned, and dereplicated output of vMAG_generation.sh
#separate out sequences into their own files for dereplicated vMAGs:
cd 14_vMAGS/all_vMAGs_QC_passed_dereplicated/vRhyme_dereplication
cat vRhyme_derep_longest_all_vMAGs_qc_passed.fa | awk '{
        if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fasta")}
        print $0 > filename
}'
mkdir vRhyme_derep_longest_all_vMAGs_qc_passed
mv *.fasta vRhyme_derep_longest_all_vMAGs_qc_passed

## RUN COVERM ON MAGS ##
coverm genome --genome-fasta-directory 10_MAPPING/MAGs_only_derep -1 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S3-paired_QC_trimmed_R_1.fastq -2 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S3-paired_QC_trimmed_R_2.fastq -x .fasta -m relative_abundance tpm rpkm mean trimmed_mean -o 10_MAPPING/coverm_MAGs_only_derep_S3_08.26.tsv
coverm genome --genome-fasta-directory 10_MAPPING/MAGs_only_derep -1 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S4-paired_QC_trimmed_R_1.fastq -2 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S4-paired_QC_trimmed_R_2.fastq -x .fasta -m relative_abundance tpm rpkm mean trimmed_mean -o 10_MAPPING/coverm_MAGs_only_derep_S4_08.26.tsv
coverm genome --genome-fasta-directory 10_MAPPING/MAGs_only_derep -1 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S5-paired_QC_trimmed_R_1.fastq -2 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S5-paired_QC_trimmed_R_2.fastq -x .fasta -m relative_abundance tpm rpkm mean trimmed_mean -o 10_MAPPING/coverm_MAGs_only_derep_S5_08.26.tsv
coverm genome --genome-fasta-directory 10_MAPPING/MAGs_only_derep -1 01_QC/*_R1.fastq.gz -2 01_QC/*_R2.fastq.gz -x .fasta -m relative_abundance tpm rpkm mean trimmed_mean -o 10_MAPPING/coverm_MAGs_only_derep_new_samples_08.26.tsv

## RUN COVERM ON vMAGS ##
coverm genome --genome-fasta-directory 14_vMAGS/all_vMAGs_QC_passed_dereplicated/vRhyme_dereplication/vRhyme_derep_longest_all_vMAGs_qc_passed -1 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S3-paired_QC_trimmed_R_1.fastq -2 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S3-paired_QC_trimmed_R_2.fastq -x .fasta -m relative_abundance tpm rpkm mean trimmed_mean -o 10_MAPPING/coverm_vMAGs_only_derep_S3_08.26.tsv
coverm genome --genome-fasta-directory 14_vMAGS/all_vMAGs_QC_passed_dereplicated/vRhyme_dereplication/vRhyme_derep_longest_all_vMAGs_qc_passed -1 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S4-paired_QC_trimmed_R_1.fastq -2 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S4-paired_QC_trimmed_R_2.fastq -x .fasta -m relative_abundance tpm rpkm mean trimmed_mean -o 10_MAPPING/coverm_vMAGs_only_derep_S4_08.26.tsv
coverm genome --genome-fasta-directory 14_vMAGS/all_vMAGs_QC_passed_dereplicated/vRhyme_dereplication/vRhyme_derep_longest_all_vMAGs_qc_passed -1 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S5-paired_QC_trimmed_R_1.fastq -2 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S5-paired_QC_trimmed_R_2.fastq -x .fasta -m relative_abundance tpm rpkm mean trimmed_mean -o 10_MAPPING/coverm_vMAGs_only_derep_S5_08.26.tsv
coverm genome --genome-fasta-directory 14_vMAGS/all_vMAGs_QC_passed_dereplicated/vRhyme_dereplication/vRhyme_derep_longest_all_vMAGs_qc_passed -1 01_QC/*_R1.fastq.gz -2 01_QC/*_R2.fastq.gz -x .fasta -m relative_abundance tpm rpkm mean trimmed_mean -o 10_MAPPING/coverm_vMAGs_only_derep_newreads_08.26.tsv

## REPEAT FOR COUNTS METRIC ##
coverm genome --genome-fasta-directory 10_MAPPING/MAGs_only_derep -1 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S3-paired_QC_trimmed_R_1.fastq -2 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S3-paired_QC_trimmed_R_2.fastq -x .fasta -m counts -o 10_MAPPING/coverm_MAGs_only_derep_counts_S3_08.26.tsv --min-covered-fraction 0
coverm genome --genome-fasta-directory 10_MAPPING/MAGs_only_derep -1 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S4-paired_QC_trimmed_R_1.fastq -2 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S4-paired_QC_trimmed_R_2.fastq -x .fasta -m counts -o 10_MAPPING/coverm_MAGs_only_derep_counts_S4_08.26.tsv --min-covered-fraction 0
coverm genome --genome-fasta-directory 10_MAPPING/MAGs_only_derep -1 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S5-paired_QC_trimmed_R_1.fastq -2 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S5-paired_QC_trimmed_R_2.fastq -x .fasta -m counts -o 10_MAPPING/coverm_MAGs_only_derep_counts_S5_08.26.tsv --min-covered-fraction 0
coverm genome --genome-fasta-directory 10_MAPPING/MAGs_only_derep -1 01_QC/*_R1.fastq.gz -2 01_QC/*_R2.fastq.gz -x .fasta -m counts -o 10_MAPPING/coverm_MAGs_only_derep_counts_new_samples_08.26.tsv --min-covered-fraction 0


coverm genome --genome-fasta-directory 14_vMAGS/all_vMAGs_QC_passed_dereplicated/vRhyme_dereplication/vRhyme_derep_longest_all_vMAGs_qc_passed -1 01_QC/*_R1.fastq.gz -2 01_QC/*_R2.fastq.gz -x .fasta -m count -o 10_MAPPING/coverm_vMAGs_only_derep_counts_newreads_08.26.tsv --min-covered-fraction 0
coverm genome --genome-fasta-directory 14_vMAGS/all_vMAGs_QC_passed_dereplicated/vRhyme_dereplication/vRhyme_derep_longest_all_vMAGs_qc_passed -1 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S3-paired_QC_trimmed_R_1.fastq -2 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S3-paired_QC_trimmed_R_2.fastq -x .fasta -m count -o 10_MAPPING/coverm_vMAGs_only_derep_counts_S3_08.26.tsv --min-covered-fraction 0
coverm genome --genome-fasta-directory 14_vMAGS/all_vMAGs_QC_passed_dereplicated/vRhyme_dereplication/vRhyme_derep_longest_all_vMAGs_qc_passed -1 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S4-paired_QC_trimmed_R_1.fastq -2 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S4-paired_QC_trimmed_R_2.fastq -x .fasta -m count -o 10_MAPPING/coverm_vMAGs_only_derep_counts_S4_08.26.tsv --min-covered-fraction 0
coverm genome --genome-fasta-directory 14_vMAGS/all_vMAGs_QC_passed_dereplicated/vRhyme_dereplication/vRhyme_derep_longest_all_vMAGs_qc_passed -1 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S5-paired_QC_trimmed_R_1.fastq -2 /workspace/mmoynihan/TR_2021/01_QC/Bloom/*S5-paired_QC_trimmed_R_2.fastq -x .fasta -m count -o 10_MAPPING/coverm_vMAGs_only_derep_counts_S5_08.26.tsv --min-covered-fraction 0
