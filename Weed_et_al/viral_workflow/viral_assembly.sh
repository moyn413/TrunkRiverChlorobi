#!/bin/bash
#assembly_megahit.sh

#make sample list
cd 00_SAMPLES
ls *_1P.gz | cut -d '-' -f 1 > sample_list

module load megahit
cd /workspace/rweed/TR_phage_redo
mkdir 02_ASSEMBLY

#indiv assemblies
for sample in $(cat sample_list); do \
     megahit --presets meta-sensitive --min-contig-len 1000 -1 01_QC/${sample}-QUALITY_PASSED_trimmed_1P.gz -2 01_QC/${sample}-QUALITY_PASSED_trimmed_2P.gz --out-prefix ${sample} -o 02_ASSEMBLY/megahit_assembly_${sample} -t 40
done

#coassembly
cd /workspace/rweed/TR_phage_redo/01_QC
cat *S4-QUALITY_PASSED_trimmed_1P.gz > concatenated_S4_1.fastq
cat *S5-QUALITY_PASSED_trimmed_1P.gz > concatenated_S5_1.fastq
cat *S4-QUALITY_PASSED_trimmed_2P.gz > concatenated_S4_2.fastq
cat *S5-QUALITY_PASSED_trimmed_2P.gz > concatenated_S5_2.fastq

cd /workspace/rweed/TR_phage_redo
megahit --presets meta-sensitive --min-contig-len 1000 -1 01_QC/concatenated_S4_1.fastq -2 01_QC/concatenated_S4_2.fastq --out-prefix coassembly_S4 -o 02_ASSEMBLY/megahit_coassembly_S4 -t 40
megahit --presets meta-sensitive --min-contig-len 1000 -1 01_QC/concatenated_S5_1.fastq -2 01_QC/concatenated_S5_2.fastq --out-prefix coassembly_S5 -o 02_ASSEMBLY/megahit_coassembly_S5 -t 40
