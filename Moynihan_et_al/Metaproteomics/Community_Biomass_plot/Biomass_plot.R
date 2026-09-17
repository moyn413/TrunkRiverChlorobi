library(plyr)
library(dplyr)
library(microViz)
library(tidyverse)
library(phyloseq)
library(RColorBrewer)

setwd("/path/")
biomass<-read.csv("Protein_Biomass", header=TRUE)

rownames(biomass)<- biomass$Accession

# Organize taxonomy and binID data
tax<-biomass %>% select(1:7)%>%  mutate_if(is.character, as.factor)

#Make abundance table
abund <- biomass %>% select(9:45)

# Okay, now I have my "otu" tables 
# And my taxonomy table, both with matching binIDs (aka otuIDs)
meta<-read.csv("MatWebBloom_Biomass_contextual_data.csv", header=TRUE)
str(meta)
rownames(meta)<- meta$ID
meta<- meta %>% 
  select(-ID)
meta<- meta %>% 
  mutate_if(is.character, as.factor)


#--------------------
#Now make phyloseq object for Biomass
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


#--------------------
# Normalize abundances
#--------------------

ps_normalize_median <- function (ps, title) {
  ps_median = median(sample_sums(ps))
  cat(sprintf("\nThe median number of reads used for normalization of %s is  %.0f", title, ps_median))
  normalize_median = function(x, t=ps_median) (if(sum(x) > 0){ t * (x / sum(x))} else {x})
  ps = transform_sample_counts(ps, normalize_median)
  cat(str_c("\nPhyloseq ",title, "\n========== \n") )
  print(ps)
}


ps.norm = ps_normalize_median(ps, "ps_norm")


=
#--------------------
# Abundance plots (Bloom B05 Samples Only) (for PAPER) 
#--------------------

#--------------------
# Make custom color palette
#--------------------
custom_palette <- c("Chlorobiales" = "#FFE100", "unbinned"= "#33322e", "Chromatiales"= "#C20088" , "Cyanobacteriales" = "#2BCE48" ,
                    "Enterobacterales"  = "#FF0010","Desulfuromonadales"= "#5EF1F2" , "Campylobacterales"= "#FFCC99" ,"Desulfobacterales"= "#F0A3FF" ,
                    "Rhodobacterales"= "#a502f0" ,"Desulfovibrionales"= "#426600" , "Chlorobiales/Chromatiales" = "#edf2a2", "Bacteroidales" = "#00998F",
                    "Thiomicrospirales" = "#4C005C", "Pseudomonadales" = "#f56042",  "Spirochaetales"  = "#993F00", "Chitinophagales"= "#c3e6e2", 
                    "Sphaerochaetales" = "#8F7C00","Kiritimatiellales" = "#990000" , "Actinomycetales"= "#eba12a", "Lachnospirales" = "#94FFB5" , "Other genera" ="lightgrey" )

#--------------------

ps.relative <- transform_sample_counts(ps.norm, function(OTU) OTU/sum(OTU))

# Make biomass abundance plot of just BloomB05 samples (primary bloom sampled during timeseries) and prebloom timepoint. Order gy sampling date. 
# Select only the 20 most abundant taxonomic Orders. 

bloom <- ps.relative %>% ps_filter(SampleType2== "Bloom05"| SampleType2=="PreBloom") %>%
  comp_barplot(
    tax_level = "Order",
    label = "Date", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    #facet_by = "SampleType2",
    tax_transform_for_plot = "identity",
    sample_order=c("TR06_PB", "TR08_B05", "TR09_B05", "TR10_B05", "TR11_B05", "TR12_B05", "TR13_B05", "TR14_B05", "TR15_B05","TR16_B05", "TR17_B05_MW"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    palette = custom_palette,
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) 


ggsave("Protein_biomass_B05_order_20_norm.pdf", plot = bloom ,
       path = "Plots_Biomass/",
       width = 10,
       height = 5)



