library(plyr)
library(dplyr)
library(microViz)
library(tidyverse)
library(phyloseq)
library(RColorBrewer)

setwd("/Users/molly/Documents/GitHub/TrunkRiverChlorobi/Moynihan_et_al/Metaproteomics/Anvio_Proteomics_Plot/GSB-TRL01_protein_abundance_by_category")

#This file has the top 250 most abundant proteins, which have been assigned broad functional categories. This assignment was done manually.
protein<-read.csv("GSB-TRL01_orgNorm_protein_functional_categories.csv", header=TRUE)

abund <- protein %>% select(3:10)
tax<-protein %>% select(2:2)

meta<-read.csv("contextual_data.csv", header=TRUE)
str(meta)
rownames(meta)<- meta$ID
meta<- meta %>% 
  select(-ID)
meta<- meta %>% 
  mutate_if(is.character, as.factor)


#--------------------
#Now make phyloseq object
abund<- as.matrix(abund)
tax<- as.matrix(tax)
class(abund)
class(tax)

#transform data to phyloseq objects
phylo_abund<- otu_table(abund, taxa_are_rows = TRUE)
phylo_TAX<- tax_table(tax)
phylo_samples<- sample_data(meta)

#and put them in one object
ps<- merge_phyloseq(phylo_abund, phylo_TAX, phylo_samples)


x1 <-ps %>% tax_transform(transformation="log10") %>% comp_barplot(
    tax_level = "Functional_category",
    label = "Date", # name an alternative variable to label axis
    n_taxa = 25, # give more taxa unique colours
    #facet_by = "SampleType2",
    tax_transform_for_plot = "identity",
    sample_order=c("TR08_B05", "TR09_B05", "TR10_B05", "TR11_B05", "TR12_B05", "TR13_B05", "TR14_B05", "TR15_B05"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
   # palette = custom_palette,
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) 


ggsave("GSB-TRL01_orgNorm_protein_category.pdf", plot = x2 ,
       path = "/path/",
       width = 10,
       height = 5)

