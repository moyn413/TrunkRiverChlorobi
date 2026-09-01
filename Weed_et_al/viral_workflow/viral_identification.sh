#!/bin/bash
#this file picks up after viral_assembly.sh
#run genomad
module load genomad
cd /workspace/rweed/TR_phage_redo
mkdir 03_CLASSIFICATION
mkdir 03_CLASSIFICATION/genomad
genomad end-to-end --cleanup -t 40 02_ASSEMBLY/megahit_coassembly_S4/coassembly_S4.contigs.fa 03_CLASSIFICATION/genomad/genomad_output_coassembly_S4 /workspace/rweed/TR_phage/03_PHAGE/genomad/genomad_db
genomad end-to-end --cleanup -t 40 02_ASSEMBLY/megahit_coassembly_S5/coassembly_S5.contigs.fa 03_CLASSIFICATION/genomad/genomad_output_coassembly_S5 /workspace/rweed/TR_phage/03_PHAGE/genomad/genomad_db

#run virsorter
module load virsorter/2.2.4-mamba
cd /workspace/rweed/TR_phage_redo
mkdir 03_CLASSIFICATION/vs2
virsorter run --keep-original-seq --seqname-suffix-off -i 02_ASSEMBLY/megahit_coassembly_S4/coassembly_S4.contigs.fa -w 03_CLASSIFICATION/vs2/vs2_output_coassembly_S4 --min-length 1000 --min-score 0.5 -j 28 all --prep-for-dramv
virsorter run --keep-original-seq --seqname-suffix-off -i 02_ASSEMBLY/megahit_coassembly_S5/coassembly_S5.contigs.fa -w 03_CLASSIFICATION/vs2/vs2_output_coassembly_S5 --min-length 1000 --min-score 0.5 -j 28 all --prep-for-dramv

#combine outputs and dereplicate
seqkit rmdup -s < <(cat genomad/genomad_output_coassembly_S5/coassembly_S5.contigs_summary/coassembly_S5.contigs_virus.fna vs2/vs2_output_coassembly_S5/final-viral-combined.fa) > coassembly_S5_viral_id_combined_fp.fa
seqkit rmdup -s < <(cat genomad/genomad_output_coassembly_S4/coassembly_S4.contigs_summary/coassembly_S4.contigs_virus.fna vs2/vs2_output_coassembly_S4/final-viral-combined.fa) > coassembly_S4_viral_id_combined_fp.fa

#run checkv
module load checkv 

inpath=/workspace/rweed/TR_phage_redo/03_CLASSIFICATION
outpath=/workspace/rweed/TR_phage_redo/03_CLASSIFICATION/CheckV
dbpath=/blastdb/checkv-db-v1.5/

#Note: needed to do some manual removing of a few duplicated sequences which resulted from trimming here (for S5 only)
checkv end_to_end $inpath/coassembly_S5_viral_id_combined_fp_edit.fa $outpath -d $dbpath -t 40
checkv end_to_end $inpath/coassembly_S4_viral_id_combined_fp.fa $outpath/S4 -d $dbpath -t 40
