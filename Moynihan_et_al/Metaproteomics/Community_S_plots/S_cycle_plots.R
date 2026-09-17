library(plyr)
library(dplyr)
library(microViz)
library(tidyverse)
library(phyloseq)
library(gridExtra)
library(RColorBrewer)


setwd("/path/")

S<-read.csv("s_cycle_all_R.csv")

attach(S)


#Make phyloseq object for plotting

# Organize taxonomy and binID data
tax<-S %>% select(Gene:Tax_bin)
binID <- S %>% select(bin_gene)
tax<- cbind(binID, tax)

rownames(tax)<- tax$bin_gene
tax<- tax %>% select(-bin_gene)

#Protein abundance data
abund<-S%>%select(TR06_PB:TR17_M01_sup)
abund<- cbind(binID, abund)

rownames(abund)<- abund$bin_gene
abund<- abund %>% select(-bin_gene)


#Load contextul data file
meta<-read.csv("contextual_data.csv", header=TRUE)
str(meta)
rownames(meta)<- meta$Sample
meta<- meta %>% 
  select(-Sample)
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
ps_gene<- merge_phyloseq(phylo_abund, phylo_TAX, phylo_samples)


#==============================================================#
# Plots for each gene, with same y-axis maximum, colored by taxa
#==============================================================#


#=======Makek Color Palette for all plots====================#
tax_df <- as.data.frame(as.matrix(tax_table(ps_gene)))  # Convert to dataframe for easier manipulation
tax_df <- tax_df[, !colnames(tax_df) %in% "Gene"]  # Remove Gene column
ps_noG<- ps_gene
tax_table(ps_noG) <- tax_table(as.matrix(tax_df))

myPal2 <- tax_palette(data = ps_noG, rank = "Order", n = 25, pal = "brewerPlus")
tax_palette_plot(myPal2)

#Make genes from the same cluster a similar color
myPal2["Oscillatoriales"] <- "#A6CEE3"
myPal2["Geitlerinematales"] <- "#1F78E4"
myPal2["Chlorobiales"] <- "#B2DF8A"
myPal2["Chromatiales"] <- "#CAB2D6"
myPal2["Enterobacterales"] <- "#FFFF99"
myPal2["Bacteroidales"] <- "#FDBF6F"
myPal2["Thiomicrospirales"] <- "#E31A1C"
myPal2["Rhodobacterales"] <- "#F1B6DA"
myPal2["Rhizobiales"] <- "#D01C8B"
myPal2["Desulfobacterales"] <- "#B15928"
myPal2["Class Cyanobacteriia"] <- "#1FF8FF"
myPal2["Desulfovibrionales"] <- "#634534"
myPal2["unbinned_web"] <- "#2b3d23"
myPal2["unbinned_mat"] <- "#2b3d24"
myPal2["unbinned_B05"] <- "#2b3d22"

tax_palette_plot(myPal2)
#=================#


# Make separate phyloseq object for main bloom (B05) samples, and microbial mat + suspended-biofilm like matrix ("web") samples
ps_B05<-ps_gene%>%ps_filter(Type != "White") %>% ps_filter(Type2 != "Bloom_B01") %>% ps_filter(Type2 != "Supernatant_B01")%>% ps_filter(Type2 != "Supernatant_B05") %>% ps_filter(Type != "Mat")%>% ps_filter(Type != "Web")
ps_MW<-ps_gene%>%ps_filter(Type != "White") %>% ps_filter(Type2 != "Bloom_B01") %>% ps_filter(Type2 != "Supernatant_B01")%>% ps_filter(Type2 != "Supernatant_B05") %>% ps_filter(Type != "Bloom")%>% ps_filter(Type != "PreBloom")%>% ps_filter(ID != "TR17_M01_DG")%>% ps_filter(ID != "TR17_M01_LG")


#==============================================================#
# Make individual plots for each sulfur protein for B05 samples, ordered by date
#==============================================================#

#----------------------------
# Sqr
#---------------------------
sqr <-ps_B05%>% subset_taxa(Gene=="sqr") %>%
  comp_barplot(
    tax_level = "Order",
    label = "ID", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    palette = myPal2,
    tax_transform_for_plot = "identity",
    sample_order=c("TR06_PB", "TR08_B05", "TR09_B05","TR10_B05","TR11_B05", "TR12_B05","TR13_B05","TR14_B05", "TR15_B05","TR16_B05","TR17_B05_MW"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) + ylim(0,0.5) + ggtitle("sqr")

#----------------------------
# fcc
#---------------------------
fcc <-ps_B05%>% subset_taxa(Gene=="fccB") %>%
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    palette = myPal2,
    tax_transform_for_plot = "identity",
    sample_order=c("TR06_PB", "TR08_B05", "TR09_B05","TR10_B05","TR11_B05", "TR12_B05","TR13_B05","TR14_B05", "TR15_B05","TR16_B05","TR17_B05_MW"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) + ylim(0,0.5) + ggtitle("fccB")


#----------------------------
# psr_phsA
#---------------------------
psr_phsA <-ps_B05%>% subset_taxa(Gene=="psr_phsA") %>%
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    palette = myPal2,
    tax_transform_for_plot = "identity",
    sample_order=c("TR06_PB", "TR08_B05", "TR09_B05","TR10_B05","TR11_B05", "TR12_B05","TR13_B05","TR14_B05", "TR15_B05","TR16_B05","TR17_B05_MW"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) + ylim(0,0.5) + ggtitle("psr_phsA")


#----------------------------
# dsrA
#---------------------------
dsrA <-ps_B05%>% subset_taxa(Gene=="dsrA") %>%
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    palette = myPal2,
    tax_transform_for_plot = "identity",
    sample_order=c("TR06_PB", "TR08_B05", "TR09_B05","TR10_B05","TR11_B05", "TR12_B05","TR13_B05","TR14_B05", "TR15_B05","TR16_B05","TR17_B05_MW"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) + ylim(0,0.5) + ggtitle("dsrA")


#----------------------------
# rdsrA
#---------------------------
rdsrA <-ps_B05%>% subset_taxa(Gene=="rdsrA") %>%
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    palette = myPal2,
    tax_transform_for_plot = "identity",
    sample_order=c("TR06_PB", "TR08_B05", "TR09_B05","TR10_B05","TR11_B05", "TR12_B05","TR13_B05","TR14_B05", "TR15_B05","TR16_B05","TR17_B05_MW"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) + ylim(0,0.5) + ggtitle("rdsrA")

#----------------------------
# aprA-oxidative
#---------------------------
aprA.ox <-ps_B05%>% subset_taxa(Gene=="aprA-ox") %>%
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    palette = myPal2,
    tax_transform_for_plot = "identity",
    sample_order=c("TR06_PB", "TR08_B05", "TR09_B05","TR10_B05","TR11_B05", "TR12_B05","TR13_B05","TR14_B05", "TR15_B05","TR16_B05","TR17_B05_MW"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) + ylim(0,0.5) + ggtitle("aprA")

#----------------------------
# aprA-reductive
#---------------------------
aprA.red <-ps_B05%>% subset_taxa(Gene=="aprA-red") %>%
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    palette = myPal2,
    tax_transform_for_plot = "identity",
    sample_order=c("TR06_PB", "TR08_B05", "TR09_B05","TR10_B05","TR11_B05", "TR12_B05","TR13_B05","TR14_B05", "TR15_B05","TR16_B05","TR17_B05_MW"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) + ylim(0,0.5) + ggtitle("aprA")


#----------------------------
# sat-ox
#---------------------------
sat.ox <-ps_B05%>% subset_taxa(Gene=="sat-ox") %>%
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    palette = myPal2,
    tax_transform_for_plot = "identity",
    sample_order=c("TR06_PB", "TR08_B05", "TR09_B05","TR10_B05","TR11_B05", "TR12_B05","TR13_B05","TR14_B05", "TR15_B05","TR16_B05","TR17_B05_MW"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) + ylim(0,0.5) + ggtitle("sat")

#----------------------------
# sat-red
#---------------------------
sat.red <-ps_B05%>% subset_taxa(Gene=="sat-red") %>%
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    palette = myPal2,
    tax_transform_for_plot = "identity",
    sample_order=c("TR06_PB", "TR08_B05", "TR09_B05","TR10_B05","TR11_B05", "TR12_B05","TR13_B05","TR14_B05", "TR15_B05","TR16_B05","TR17_B05_MW"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) + ylim(0,0.5) + ggtitle("sat")


#----------------------------
# sdo
#---------------------------
sdo <-ps_B05%>% subset_taxa(Gene=="sdo") %>%
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    palette = myPal2,
    tax_transform_for_plot = "identity",
    sample_order=c("TR06_PB", "TR08_B05", "TR09_B05","TR10_B05","TR11_B05", "TR12_B05","TR13_B05","TR14_B05", "TR15_B05","TR16_B05","TR17_B05_MW"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) + ylim(0,0.5) + ggtitle("sdo")


#----------------------------
# soxA
#---------------------------
soxA <-ps_B05%>% subset_taxa(Gene=="soxA") %>%
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    palette = myPal2,
    tax_transform_for_plot = "identity",
    sample_order=c("TR06_PB", "TR08_B05", "TR09_B05","TR10_B05","TR11_B05", "TR12_B05","TR13_B05","TR14_B05", "TR15_B05","TR16_B05","TR17_B05_MW"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) + ylim(0,0.5) + ggtitle("soxA")


x.order <-grid.arrange(sqr, fcc, psr_phsA, sdo, soxA, rdsrA, 
                       dsrA, sat.ox, sat.red, aprA.red, aprA.ox , nrow = 1)

ggsave("Scycle_bygene_B05.pdf", plot = x.order, 
       path = "/path/", limitsize = FALSE,
       width = 55,
       height = 5)




#==============================================================#
# Make individual plots for each sulfur protein for Mat and Suspended Biofilm (" web" ) samples
#==============================================================#

#----------------------------
# Sqr
#---------------------------
sqr <-ps_MW%>% subset_taxa(Gene=="sqr") %>%
  comp_barplot(
    tax_level = "Order",
    label = "ID", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    palette = myPal2,
    tax_transform_for_plot = "identity",
    sample_order=c("TR12_W_P", "TR13_W_P", "TR12_M01","TR13_M01","TR14_M01", "TR16_M01","TR17_M01"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) + ylim(0,0.5) + ggtitle("sqr")

#----------------------------
# fcc
#---------------------------
fcc <-ps_MW%>% subset_taxa(Gene=="fccB") %>%
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    palette = myPal2,
    tax_transform_for_plot = "identity",
    sample_order=c("TR12_W_P", "TR13_W_P", "TR12_M01","TR13_M01","TR14_M01", "TR16_M01","TR17_M01"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) + ylim(0,0.5) + ggtitle("fccB")


#----------------------------
# psr_phsA
#---------------------------
psr_phsA <-ps_MW%>% subset_taxa(Gene=="psr_phsA") %>%
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    palette = myPal2,
    tax_transform_for_plot = "identity",
    sample_order=c("TR12_W_P", "TR13_W_P", "TR12_M01","TR13_M01","TR14_M01", "TR16_M01","TR17_M01"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) + ylim(0,0.5) + ggtitle("psr_phsA")


#----------------------------
# dsrA
#---------------------------
dsrA <-ps_MW%>% subset_taxa(Gene=="dsrA") %>%
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    palette = myPal2,
    tax_transform_for_plot = "identity",
    sample_order=c("TR12_W_P", "TR13_W_P", "TR12_M01","TR13_M01","TR14_M01", "TR16_M01","TR17_M01"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) + ylim(0,0.5) + ggtitle("dsrA")


#----------------------------
# rdsrA
#---------------------------
rdsrA <-ps_MW%>% subset_taxa(Gene=="rdsrA") %>%
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    palette = myPal2,
    tax_transform_for_plot = "identity",
    sample_order=c("TR12_W_P", "TR13_W_P", "TR12_M01","TR13_M01","TR14_M01", "TR16_M01","TR17_M01"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) + ylim(0,0.5) + ggtitle("rdsrA")

#----------------------------
# aprA-oxidative
#---------------------------
aprA.ox <-ps_MW%>% subset_taxa(Gene=="aprA-ox") %>%
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    palette = myPal2,
    tax_transform_for_plot = "identity",
    sample_order=c("TR12_W_P", "TR13_W_P", "TR12_M01","TR13_M01","TR14_M01", "TR16_M01","TR17_M01"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) + ylim(0,0.5) + ggtitle("aprA")

#----------------------------
# aprA-reductive
#---------------------------
aprA.red <-ps_MW%>% subset_taxa(Gene=="aprA-red") %>%
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    palette = myPal2,
    tax_transform_for_plot = "identity",
    sample_order=c("TR12_W_P", "TR13_W_P", "TR12_M01","TR13_M01","TR14_M01", "TR16_M01","TR17_M01"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) + ylim(0,0.5) + ggtitle("aprA")


#----------------------------
# sat-ox
#---------------------------
sat.ox <-ps_MW%>% subset_taxa(Gene=="sat-ox") %>%
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    palette = myPal2,
    tax_transform_for_plot = "identity",
    sample_order=c("TR12_W_P", "TR13_W_P", "TR12_M01","TR13_M01","TR14_M01", "TR16_M01","TR17_M01"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) + ylim(0,0.5) + ggtitle("sat")

#----------------------------
# sat-red
#---------------------------
sat.red <-ps_MW%>% subset_taxa(Gene=="sat-red") %>%
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    palette = myPal2,
    tax_transform_for_plot = "identity",
    sample_order=c("TR12_W_P", "TR13_W_P", "TR12_M01","TR13_M01","TR14_M01", "TR16_M01","TR17_M01"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) + ylim(0,0.5) + ggtitle("sat")


#----------------------------
# sdo
#---------------------------
sdo <-ps_MW%>% subset_taxa(Gene=="sdo") %>%
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    palette = myPal2,
    tax_transform_for_plot = "identity",
    sample_order=c("TR12_W_P", "TR13_W_P", "TR12_M01","TR13_M01","TR14_M01", "TR16_M01","TR17_M01"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) + ylim(0,0.5) + ggtitle("sdo")


#----------------------------
# soxA
#---------------------------
soxA <-ps_MW%>% subset_taxa(Gene=="soxA") %>%
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    palette = myPal2,
    tax_transform_for_plot = "identity",
    sample_order=c("TR12_W_P", "TR13_W_P", "TR12_M01","TR13_M01","TR14_M01", "TR16_M01","TR17_M01"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) + ylim(0,0.5) + ggtitle("soxA")


x.order <-grid.arrange(sqr, fcc, psr_phsA, sdo, soxA, rdsrA, 
                       dsrA, sat.ox, sat.red, aprA.red, aprA.ox , nrow = 1)

ggsave("Scycle_bygene_MW.pdf", plot = x.order, 
       path = "/path/", limitsize = FALSE,
       width = 55,
       height = 5)








