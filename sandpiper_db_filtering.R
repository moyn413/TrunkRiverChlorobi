
#combine sandpiper1.1.0.globdb.with_extras_species_50p with sandpiper1.1.0.per_acc_summary.csv
per_acc_summary <- read.csv("sandpiper1.1.0.per_acc_summary.csv", header = TRUE)
species_50p <- read.csv("sandpiper1.1.0.globdb.with_extras_species_50p.tsv", header = FALSE, sep = "\t")
species_50p_metadata <- merge(species_50p, per_acc_summary, by.x = "V1", by.y = "sample")
#cleanup:
species_50p_metadata <- species_50p_metadata[-5]
library(dplyr)
species_50p_metadata <- species_50p_metadata %>% rename(sample = V1, coverage1 = V2, coverage2 = V3, relative_abundance = V4, taxonomic_classification = V6)

#remove anything host associated:
species_50p_metadata_environmental <- species_50p_metadata %>% filter(host_or_not != "host")
#there are still a lot of things that are predicted as hosts, so also filter by the predicted column
species_50p_metadata_environmental <- species_50p_metadata_environmental %>% filter(prediction != "host")

#only take the environments that are potentially *aquatic blooms* of any kind.
#I would summarize my criteria by saying that I selected natural (not associated with wastewater processing or fermentation or any industrial process, non-host affiliated, aquatic samples (excluding sedimentary samples). I also excluded hot springs since they have been shown to show interesting defense patterns which likely have to do with the high heat. 
aquatic_environments <- c("lake water metagenome", "brine metagenome", "aquatic metagenome", "riverine metagenome", "groundwater metagenome", "salt lake metagenome", "pond metagenome", "estuary metagenome", "soda lake metagenome", "marsh metagenome", "marine metagenome", "freshwater metagenome", "aquifer metagenome", "mangrove metagenome", "hypersaline lake metagenome", "salt marsh metagenome", "seawater metagenome", "wetland metagenome", "marine plankton metagenome")
species_50p_metadata_aquatic <- species_50p_metadata_environmental %>% filter(organism %in% aquatic_environments)
#1599 samples, only 776 unique taxa.

#export and check manually:
write.table(species_50p_metadata_aquatic, file='species_50p_metadata_aquatic.tsv', quote=FALSE, sep='\t', row.names=FALSE)

#complete list of filtering (some of which was done via manual check)
#1. Has greater than 50% relative abundance
#2. Is in one of the following categories: c("lake water metagenome", "brine metagenome", "aquatic metagenome", "riverine metagenome", "groundwater metagenome", "salt lake metagenome", "pond metagenome", "estuary metagenome", "soda lake metagenome", "marsh metagenome", "marine metagenome", "freshwater metagenome", "aquifer metagenome", "mangrove metagenome", "hypersaline lake metagenome", "salt marsh metagenome", "seawater metagenome", "wetland metagenome", "marine plankton metagenome")
#3. Has a prokaryote fraction greater than 50%
#4. Everything that passes these criteria is checked manually to make sure it is a true environmental sample, not from a hot spring, and not host associated (except a few pathogens which are also abundant in the environment)
#5. Everything is also checked to make sure that it has high relative abundance in more than one single instance
#from these criteria, I recovered only 58 samples total
