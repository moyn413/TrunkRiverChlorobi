#!/bin/bash

##GET LIST OF GLOBDBIDS##
wget https://fileshare.lisc.univie.ac.at/globdb/globdb_r226/globdb_r226_taxonomy.tsv.gz

cut -f 1 bloom_former_review.tsv > bloomers.txt

while IFS= read -r bloomer; do
    result=$(grep -F "$bloomer" globdb_r226_taxonomy.tsv | paste -sd ';' -)
    [[ -z $result ]] && result="NA"
    printf "%s\t%s\n" "$bloomer" "$result"
done < bloomers.txt > glob_ids_bloomers.tsv

#ok so many of these have NAs, and some have returned multiple things. 
#manually curate and save in Bloom_former_review_final.txt


##GET BLOOMER GENOME LOCATIONS##
cut -f 11 Bloom_former_review_final.txt > glob_ids.txt
grep -f glob_ids.txt /usr/local/tmp/globdb_r226_genome_fasta.tar.gz-catalog.txt | tee globdb_locations.txt

tar --verbose --extract --files-from globdb_locations.txt \
-f /usr/local/tmp/globdb_r226_genome_fasta.tar.gz -C bloomers_genomes
#this saves them all within the file structure unfortunately
for file in $(cat globdb_locations.txt); do mv $file bloomers_genomes; done

##RUN DEFENSE FINDER##
module load miniconda
source activate /workspace/rweed/TR_phage_redo/defensefinder_env/

#run:
for genome in $(ls bloomers_genomes/*.fa.gz); do \
	defense-finder run $genome -o df_bloomers_genomes
done

#combine all of these in output file:
awk '
    FNR == 1 {
        bin = FILENAME
        gsub(/.*\/|_defense_finder_systems\.tsv$/,"",bin)
        next
    }
    { print bin "\t" $0 }
' *_defense_finder_systems.tsv > combined_defense_finder_systems.tsv

#count the number of defense systems
for id in $(cat glob_ids.txt); do \
	grep $id -c combined_defense_finder_systems.tsv >> number_of_systems
done

paste glob_ids.txt number_of_systems > systems_by_id

##RANDOM SUBSETS OF NON-BLOOMERS##
#get the catalog excluding the bloom former IDs
diff -13 globdb_locations.txt /usr/local/tmp/globdb_r226_genome_fasta.tar.gz-catalog.txt > globdb_r226_genome_fasta.tar.gz-catalog_no_bloomers.txt

for i in {1..19}; do \
	shuf -n 20 globdb_r226_genome_fasta.tar.gz-catalog_no_bloomers.txt | tee temp
	cut -c3- temp > rand_subset_locations
	rm temp
	tar --verbose --extract --files-from rand_subset_locations -f /usr/local/tmp/globdb_r226_genome_fasta.tar.gz
	#this saves them all within the file structure unfortunately
	mkdir genomes
	for file in $(cat rand_subset_locations); do mv $file genomes; done
	rm -r globdb_r226_genome_fasta
	#run defense finder on all files:
	module load miniconda
	source activate /workspace/rweed/TR_phage_redo/defensefinder_env/
	#run:
	for genome in $(ls genomes/*.fa.gz); do \
		defense-finder run $genome -o df_genomes
	done
	
	#summarize:
	cd df_genomes
	awk '
	    FNR == 1 {
	        bin = FILENAME
	        gsub(/.*\/|_defense_finder_systems\.tsv$/,"",bin)
	        next
	    }
	    { print bin "\t" $0 }
	' *_defense_finder_systems.tsv > combined_defense_finder_systems.tsv
	
	mv combined_defense_finder_systems.tsv ..
	cd ..
	#get just IDs:
	cut -d / -f 3 rand_subset_locations > files
	cut -d . -f 1 files > glob_ids
	
	#count the number of defense systems
	for id in $(cat glob_ids); do \
		grep $id -c combined_defense_finder_systems.tsv >> number_of_systems; \
	done
	
	paste glob_ids number_of_systems > systems_by_id_${i}
	
	rm glob_ids number_of_systems combined_defense_finder_systems.tsv files rand_subset_locations
	rm -r genomes df_genomes
done

##GET GENOME STATS##
#get genome stats file from globdb in order to normalize by length
wget https://fileshare.lisc.univie.ac.at/globdb/globdb_r226/globdb_r226_tax_plus_stats.tsv.gz