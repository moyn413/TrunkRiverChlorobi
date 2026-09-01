#!/bin/bash

#get list of putative bloom formers from sandpiper:
singlem summarise --input-taxonomic-profile sandpiper1.1.0.globdb.csv --output-taxonomic-profile-with-extras sandpiper1.1.0.globdb.with_extras.tsv --metapackage GlobDB_r226.metapackage_v1.smpkg

#filter database:
#has species level ID
awk '$5=="species" {print $0}' sandpiper1.1.0.globdb.with_extras.tsv > sandpiper1.1.0.globdb.with_extras_species.tsv
#57213646
#rel abund over 50%
awk '$4>= 50 {print $0}' sandpiper1.1.0.globdb.with_extras_species.tsv > sandpiper1.1.0.globdb.with_extras_species_50p.tsv
#79979
