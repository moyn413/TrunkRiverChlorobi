#!/bin/bash
#this file picks up from viral_identification.sh
#Note: S#_QC_passed_contigs_* files were generated manually from the CheckV results in Excel using the filtering standards outlined in the text.
#S4
#make files for subsetting:
S4_QC_passed_contigs_10kb_min.txt #listed as not containing a prophage, 171 sequences
S4_QC_passed_contigs_10kb_min_proviral.txt #listed as containing a prophage, 16 sequences

#now subset from original sequence file:
seqtk subseq ../../coassembly_S4_viral_id_combined_fp.fa S4_QC_passed_contigs_10kb_min.txt > S4_QC_passed_contigs_10kb_min.fasta
seqtk subseq ../../coassembly_S4_viral_id_combined_fp.fa S4_QC_passed_contigs_10kb_min_proviral.txt > S4_QC_passed_contigs_10kb_min_proviral.fasta

#now 5kb min:
#make files for subsetting:
S4_QC_passed_contigs_5kb_min.txt #listed as not containing a prophage, 498 sequences
S4_QC_passed_contigs_5kb_min_proviral.txt #listed as containing a prophage, 22 sequences

seqtk subseq ../../coassembly_S4_viral_id_combined_fp.fa S4_QC_passed_contigs_5kb_min.txt > S4_QC_passed_contigs_5kb_min.fasta
seqtk subseq ../../coassembly_S4_viral_id_combined_fp.fa S4_QC_passed_contigs_5kb_min_proviral.txt > S4_QC_passed_contigs_5kb_min_proviral.fasta

#S5
#make files for subsetting:
S5_QC_passed_contigs_10kb_min.txt #listed as not containing a prophage, 2916 sequences
S5_QC_passed_contigs_10kb_min_proviral.txt #listed as containing a prophage, 71 sequences
S5_QC_passed_contigs_5kb_min.txt #listed as not containing a prophage, 9769 sequences
S5_QC_passed_contigs_5kb_min_proviral.txt #listed as containing a prophage, 247 sequences

#now subset from original sequence file:
seqtk subseq ../../coassembly_S5_viral_id_combined_fp_edit.fa S5_QC_passed_contigs_10kb_min.txt > S5_QC_passed_contigs_10kb_min.fasta
seqtk subseq ../../coassembly_S5_viral_id_combined_fp_edit.fa S5_QC_passed_contigs_10kb_min_proviral.txt > S5_QC_passed_contigs_10kb_min_proviral.fasta
seqtk subseq ../../coassembly_S5_viral_id_combined_fp_edit.fa S5_QC_passed_contigs_5kb_min.txt > S5_QC_passed_contigs_5kb_min.fasta
seqtk subseq ../../coassembly_S5_viral_id_combined_fp_edit.fa S5_QC_passed_contigs_5kb_min_proviral.txt > S5_QC_passed_contigs_5kb_min_proviral.fasta


#then combine co-assemblies
cat S5_QC_passed_contigs_10kb_min.fasta S5_QC_passed_contigs_10kb_min_proviral.fasta > S5_QC_passed_contigs_10kb_min_all.fasta
cat S5_QC_passed_contigs_5kb_min.fasta S5_QC_passed_contigs_5kb_min_proviral.fasta > S5_QC_passed_contigs_5kb_min_all.fasta

cat S4_QC_passed_contigs_10kb_min.fasta S4_QC_passed_contigs_10kb_min_proviral.fasta > S4_QC_passed_contigs_10kb_min_all.fasta
cat S4_QC_passed_contigs_5kb_min.fasta S4_QC_passed_contigs_5kb_min_proviral.fasta > S4_QC_passed_contigs_5kb_min_all.fasta

#combine size fractions:
cat S5/S5_QC_passed_contigs_10kb_min_all.fasta S4/S4_QC_passed_contigs_10kb_min_all.fasta > all_QC_contigs_10kb_min.fasta
cat S5/S5_QC_passed_contigs_5kb_min_all.fasta S4/S4_QC_passed_contigs_5kb_min_all.fasta > all_QC_contigs_5kb_min.fasta


