
conda activate anvio-dev


# extract gene calls and annotations from gene bank file from Geneious Prime (where I have protein names and annotations of the final refined, aligned MAG)
anvi-script-process-genbank -i Prosthecochloris_GSB-TRL01_refined_aligned.gb \
                            --output-gene-calls Prosthecochloris_GSB-TRL01_genes.tsv \
                            --output-functions Prosthecochloris_GSB-TRL01 _functions.tsv \
                            --output-fasta Prosthecochloris_GSB-TRL01.fa \
                            --annotation-source prodigal --include-locus-tags-as-functions


# make contigs db with aligned genome file and import gene calls and genes from gene bank files generated above
anvi-gen-contigs-database -f Prosthecochloris_GSB-TRL01.fa \
                          -o contigs.db \
                          --external-gene-calls \
                          Prosthecochloris_GSB-TRL01_genes.tsv \
                          --split-length -1

                            
# import functions that were extracted from gene bank file from Geneious Prime. Include hypothetical proteins.   
anvi-import-functions -c contigs.db \
                      -i Prosthecochloris_GSB-TRL01_functions.tsv


# run default anvi'o HMMs
anvi-run-hmms -c contigs.db  --num-threads 8

# run NCBI COGs
anvi-run-ncbi-cogs -c contigs.db   --num-threads 8

# run kegg K0fams
anvi-run-kegg-kofams -c contigs.db -T 8

# run pfams
anvi-run-pfams -c contigs.db -T 8



# Visualize aligned MAG with its proteome in anvi-interactive 
# Using file that has broad functional categories for the top 230 proteins. 
# These were assigned manually based on annotations and can be found in columns ProteinCategory1-3.
# Also includes hypothetical proteins. 

anvi-interactive --manual-mode -d Prosthecochloris_GSB-TRL01 _functional_categories.txt -p PROFILE.db




