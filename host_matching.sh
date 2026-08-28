#!/bin/bash
#this file uses the output of vMAG_generation.sh and {Molly's binning file}

##ADD BINS TO DATABASE##
#run CheckM
checkm lineage_wf 06_HOSTS/all_TR_bins 06_HOSTS/all_TR_bins

cd 06_HOSTS/CheckM_all_TR
awk '$2 >= 50 && $3 <= 10' quality_report.tsv > medium_quality_and_up.tsv
#299 bins passed these quality metrics

cd ../..
mkdir 06_HOSTS/medium_quality_and_up_bins
cp 11_NANOPORE/polished_contigs_separated/contig1_20_concat.fasta 06_HOSTS/medium_quality_and_up_bins/contig1_20_concat.fa

cd 06_HOSTS
cut -f1 medium_quality_and_up.tsv > medium_and_up_list

list=`cat medium_and_up_list`
for bin in $list; do cp all_TR_bins/${bin}.fa medium_quality_and_up_bins; done

#run gtdbtk
module load gtdbtk
gtdbtk de_novo_wf --genome_dir  medium_quality_and_up_bins/ --bacteria --outgroup_taxon p__Chloroflexota --out_dir  medium_quality_and_up_bins_GTDB-tk_results_02.26/ --cpus 48 --force --extension fa
gtdbtk de_novo_wf --genome_dir  medium_quality_and_up_bins/ --archaea --outgroup_taxon p__Undinarchaeota --out_dir  medium_quality_and_up_bins_GTDB-tk_results_02.26/ --cpus 48 --force --extension fa

#add bins to iphop db:
module load iphop
clusterize -n 48 -m lweed@uchicago.edu -log iphop_add_to_db_02.26.retry.log iphop add_to_db --fna_dir /workspace/rweed/TR_phage_redo/06_HOSTS/medium_quality_and_up_bins --gtdb_dir /workspace/rweed/TR_phage_redo/06_HOSTS/medium_quality_and_up_bins_GTDB-tk_results_02.26/ --out_dir workspace/rweed/TR_phage_redo/UTILS/Feb_2026_w_TR_hosts --db_dir /blastdb/iphop-db-aug23-rw/Aug_2023_pub_rw

##RUN IPHOP##
#iphop doesn't like my longer NNN linkers, so shorten them
link_bin_sequences.py -i 14_vMAGS/vMAGS_binned_and_unbinned/vRhyme_bins -o 14_vMAGS/linked_bins_100N -n 100
#combine into one fasta file:
cd linked_bins_100N
cat vRhyme* >> all_qc_derep_vMAGs_100N.fasta
#add the rest of the qc viral contigs:
mv all_qc_derep_vMAGs_100N.fasta all_vMAGs_QC_passed_dereplicated/vRhyme_dereplication/vRhyme_derep_longest_all_vMAGs_qc_passed
cd all_vMAGs_QC_passed_dereplicated/vRhyme_dereplication/vRhyme_derep_longest_all_vMAGs_qc_passed
cat k141* >> all_qc_derep_vMAGs_100N.fasta
mv all_qc_derep_vMAGs_100N.fasta vRhyme_derep_longest_all_vMAGs_qc_passed_100N_linked.fa

#run iphop
clusterize -n 48 -m lweed@uchicago.edu -log iphop_predict_02.2026.log iphop predict --fa_file /workspace/rweed/TR_phage_redo/14_vMAGS/all_vMAGs_QC_passed_dereplicated/vRhyme_dereplication/vRhyme_derep_longest_all_vMAGs_qc_passed_100N_linked.fa --db_dir /workspace/rweed/TR_phage_redo/UTILS/Feb_2026_w_TR_hosts --out_dir /workspace/rweed/TR_phage_redo/06_HOSTS/iphop_all_vMAGs_qc_passed_100N_linked_02.2026

