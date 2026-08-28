#!/bin/bash

#moving forward I am using the 5kb minimum file only
#this file picks up from viral_QC_and_combine.sh

##DEREPLICATE##
#remove duplicates (note: this step is probably not necessary, but it is what I did)
seqkit rmdup -s < all_QC_contigs_5kb_min.fasta > all_QC_contigs_5kb_min_rmdup.fasta
#dereplicate using sequence identity
cd-hit-est -i all_QC_contigs_5kb_min_rmdup.fasta -o all_QC_contigs_5kb_min_rmdup_derep.fasta -c 0.95 -aS 0.85

##PREP FOR BINNING##
#rerun CheckV (this step is not necessary, but it is what I did; note running on rmdup file rather than derep file)
checkv end_to_end all_QC_contigs_5kb_min_rmdup.fasta CheckV_redo_5kb -d /blastdb/checkv-db-v1.5/ -t 40

#exclude complete genomes and prophages from viral binning
#get complete genomes
cut -f1 /workspace/rweed/TR_phage_redo/03_CLASSIFICATION/CheckV/CheckV_redo_5kb/complete_genomes.tsv > excluded_list

#get prophages
cd /workspace/rweed/TR_phage_redo/03_CLASSIFICATION/CheckV/CheckV_redo_5kb/
grep "^>" proviruses.fna > prophage_headers
sed 's/_[^_]*$//' prophage_headers > prophage_headers_2
#breaking down the sed: 
# _ matches underscore, [^_]* matches zero or more characters that are not underscores ([^_] is the negation of an underscore), $ means that the position of the pattern to be matched is at the end of the line
#now remove > at the beginning of the line:
sed 's/>//' prophage_headers_2 > prophage_headers_3

cat prophage_headers_3 excluded_list > excluded_list_final

#remove extra info from the header 
reformat.sh in=all_QC_contigs_5kb_min_rmdup_derep.fasta  out=all_QC_contigs_5kb_min_rmdup_derep_trim.fasta  trd=t
#now extract only sequences not in list
python /workspace/rweed/Beacon/Beacon_2023_HF_MG/UTILS/faSomeRecords.py -f all_QC_contigs_5kb_min_rmdup_derep_trim.fasta --list CheckV_redo_5kb/excluded_list_final -o all_QC_contigs_5kb_min_for_vRhyme.fasta --exclude

##BIN##
#need to rename read files so that they match what vRhyme is expecting:
cd 01_QC
for file in ./*P.gz ; do mv "$file" "${file//P.gz/.fastq.gz}" ; done

vRhyme -i 03_CLASSIFICATION/CheckV/all_QC_contigs_5kb_min_for_vRhyme.fasta -r 01_QC/*.fastq.gz -t 48 -o 14_vMAGS/vRhyme

#link bins with NNNs for downstream analyses:
cd 14_vMAGS/vRhyme
link_bin_sequences.py -i vRhyme_best_bins_fasta -o linked_bins

#recombine with complete contigs, prophages, and unbinned contigs
#remove binned sequences from original input fasta (not the one that has prophages and complete genomes removed)
cut -f1 vRhyme_best_bins.0.membership.tsv > binned_contigs
python /workspace/rweed/Beacon/Beacon_2023_HF_MG/UTILS/faSomeRecords.py -f all_QC_contigs_5kb_min_rmdup_derep_trim.fasta --list ../../14_vMAGS/vRhyme/binned_contigs -o all_QC_contigs_5kb_min_unbinned.fasta --exclude

#combine all linked bins into one fasta file:
cat *.fasta >> all_vMAGs.fasta
cat vRhyme/linked_bins/all_vMAGs.fasta ../03_CLASSIFICATION/CheckV/all_QC_contigs_5kb_min_unbinned.fasta > vMAGs_unbinned_combined.fasta

##QUALITY FILTER##
#run checkV again
checkv end_to_end vMAGs_unbinned_combined.fasta

#take only contigs longer than 10kb or categorized as high quality or complete
awk '$2 >= 10000 || $8 == "High-quality" || $8 == "Complete"'  quality_summary.tsv > qc_passed_qs.tsv

##DEREPLICATE##
#subset vMAGs_unbinned_combined.fasta to take only QC passed contigs
seqtk subseq vMAGs_unbinned_combined.fasta qc_passed_list > all_vMAGs_qc_passed.fasta

#use vRhyme to dereplicate
vRhyme -i all_vMAGs_qc_passed.fasta -t 24 -o all_vMAGs_QC_passed_dereplicated/ --derep_only --method longest

##ANNOTATE AND GET TAX##
#annotate
module load cenote-taker3/
cenotetaker3 -c /workspace/rweed/TR_phage_redo/14_vMAGS/all_vMAGs_QC_passed_dereplicated/vRhyme_dereplication/vRhyme_derep_longest_all_vMAGs_qc_passed.fa -r ct3_all_vmags -p F -am T -t 40

#re-run genomad for taxonomy
genomad annotate --cleanup -t 40 /workspace/rweed/TR_phage_redo/14_vMAGS/all_vMAGs_QC_passed_dereplicated/vRhyme_dereplication/vRhyme_derep_longest_all_vMAGs_qc_passed.fa  /workspace/rweed/TR_phage_redo/14_vMAGS/all_vMAGs_QC_passed_genomad_annotate /workspace/rweed/TR_phage_redo/UTILS/genomad_db > genomad_annotate_vmags.log
