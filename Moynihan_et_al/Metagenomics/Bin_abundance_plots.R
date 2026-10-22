library(plyr)
library(dplyr)
library(microViz)
library(tidyverse)
library(phyloseq)
library(RColorBrewer)


setwd("~/Metagenomics")

cover<-read.csv("coverm_genomes.csv", header=TRUE)

#Remove unmapped row (deal with relative abundance later if needed)
cover<- cover %>% na.omit()

colnames(cover)[1]<- "binID"
rownames(cover)<- cover$binID

#Separate by size by metric
# rel.abund <- cover %>%  dplyr::select(contains("Abundance", ignore.case = TRUE)) # Ignore this one for now because of the high unmapped percentage 
rpkm <- cover %>%  dplyr::select(contains("RPKM", ignore.case = TRUE))
trim.mean <- cover %>%  dplyr::select(contains("Trimmed", ignore.case = TRUE))
tpm <- cover %>%  dplyr::select(contains("TPM", ignore.case = TRUE))
abund1 <- cover %>%  dplyr::select(contains("Relative", ignore.case = TRUE))

# #Export abund1 as csv file
# write.csv(abund1, "/Users/molly/Dropbox/POSTDOC_MBL/Data/Metagenomes/TR2021_Metagenomes/TR_2021_B05/coverm/relabundance.csv")

#Standardize column sample names
names(rpkm) <- sub(".RPKM", "", names(rpkm))
names(trim.mean) <- sub(".Trimmed.Mean", "", names(trim.mean))
names(tpm) <- sub(".TPM", "", names(tpm))
names(abund1) <- sub(".Relative.Abundance", "", names(abund1))

# Organize taxonomy and binID data
tax<-cover %>% select(2:9)%>%  mutate_if(is.character, as.factor)
binID <- cover %>% select(1:1)
tax<- cbind(binID, tax)

rownames(tax)<- tax$binID
tax<- tax %>% select(-binID)

# Okay, now I have my "otu" tables in the form of rpkm, trimmed mean, tpm values, and adjusted abundance
# And my taxonomy table, both with matching binIDs (aka otuIDs)
meta<-read.csv("contextual_data.csv", header=TRUE)
str(meta)
rownames(meta)<- meta$Sample
meta<- meta %>% 
  select(-Sample)
meta<- meta %>% 
  mutate_if(is.character, as.factor)


#make proper levels of the factor so that they are in the wanted order in figures
meta$Date<- factor(meta$Date, levels = c("26/8/21", "29/8/21", "30/8/21", "31/8/21", "1/9/21", "21/10/21"))
meta$Date2<- factor(meta$Date2, levels = c("26-Aug-21", "29-Aug-21", "30-Aug-21", "31-Aug-21", "1-Sep-21", "21-Oct-21"))
meta$SamplingDay<- factor(meta$SamplingDay, levels = c("TR08", "TR11", "TR12", "TR13", "TR14", "TR17"))
meta$Size<- factor(meta$Size, levels = c("S1", "S2", "S3", "S1-S3","S4" ,"S5"))
meta$Size2<- factor(meta$Size2, levels = c("Above 5um", "2-5um", "0.2-2um", "Above 0.2um", "0.1-0.2um", "0.025-0.1um"))


#--------------------
#Now make phyloseq object for RPKM
rpkm<- as.matrix(rpkm)
tax<- as.matrix(tax)
class(rpkm)
class(tax)


#transform data to phyloseq objects
phylo_RPKM<- otu_table(rpkm, taxa_are_rows = TRUE)
phylo_TAX<- tax_table(tax)
phylo_samples<- sample_data(meta)

#and put them in one object
phylo_rkpm<- merge_phyloseq(phylo_RPKM, phylo_TAX, phylo_samples)

#--------------------
#Now make phyloseq object for just relative abundance
abund1<- as.matrix(abund1)
tax<- as.matrix(tax)
class(abund1)
class(tax)


#transform data to phyloseq objects
phylo_abund1<- otu_table(abund1, taxa_are_rows = TRUE)
phylo_TAX<- tax_table(tax)
phylo_samples<- sample_data(meta)

#and put them in one object
phylo_abund1<- merge_phyloseq(phylo_abund1, phylo_TAX, phylo_samples)

#--------------------
# Abundance plots for main figure paper (coverM based)
#--------------------

custom_palette <- c("Chlorobiales_semibin_bin.163" = "#FFE100", 
                    "Xanthomonadales_semibin_bin.302" = "#a57df0",
                    "Chromatiales_semibin_bin.271"= "#C20088",
                    "Pseudomonadales_semibin_bin.301"  = "#f56042",
                    "Desulfuromonadales_semibin_bin.206"   = "#5EF1F2" ,
                    "Enterobacterales_A_maxbin_bin.10"  = "#FF0010",
                    "Desulfuromonadales_semibin_bin.87"      = "#0c64e8" ,
                    "Desulfuromonadales_maxbin_bin.48_sub"      = "#588ee0" ,
                    "Enterobacterales_maxbin_bin.169"        = "#fa4c07" ,
                    "Bacteroidales_metabat_bin.124"   =  "#00998F",
                    "Desulfovibrionales_semibin_bin.137_sub"  = "#426600" ,
                    "Desulfuromonadales_semibin_bin.220" = "#37d4f0",
                    "Exiguobacterales_maxbin_bin.218" = "#e4ffa1",
                    "Rhodobacterales_semibin_bin.246" = "#a502f0" ,
                    "Sphaerochaetales_maxbin_bin.27_sub" =  "#8F7C00",
                    "Peptostreptococcales_metabat_bin.181"= "#808080" ,
                    "Peptostreptococcales_metabat_bin.45" = "#dbd4d3",
                    "Moranbacterales_metabat_bin.135" ="#354502" ,
                    "Bacteroidales_metabat_bin.90" = "#afc9c8", 
                    "Chromatiales_semibin_bin.64"= "#f250c2",
                    "other" ="lightgrey" )

# Selecting S4 size range (0.1-0.2um filter fraction)
bin.abund.S4= phylo_abund1  %>% ps_filter(Size == "S4") %>% ps_filter(SamplingDay != "TR17") %>%
  comp_barplot(
    tax_level = "Order_Bin",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    facet_by = "Size",
    tax_transform_for_plot = "identity",
    sample_order=c( "TR08_B05_S4",  "TR11_B05_S4", "TR12_B05_S4", "TR13_B05_S4","TR14_B05_S4"), 
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    palette = custom_palette,
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) 


# ggsave("Abundance_S4_bins_20.pdf", plot = x2 ,
#        path = "plots/",
#        width = 10,
#        height = 5)


#--------------------
# Plots for Supplement from all filter fractions
#--------------------

# S1,  5um filter
bin.abund.S1= phylo_abund1  %>% ps_filter(Size == "S1") %>% 
  comp_barplot(
    tax_level = "Order_Bin",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    facet_by = "Size",
    tax_transform_for_plot = "identity",
    sample_order=c( "TR08_B05_S1",  "TR11_B05_S1", "TR14_B05_S1"), 
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    #palette = custom_palette,
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) 


# ggsave("Abundance_S1_bins_20.pdf", plot = bin.abund.S1 ,
#        path = "plots/",
#        width = 10,
#        height = 5)


#--------------------

# S2,  2-5um filter fraction
bin.abund.S2= phylo_abund1  %>% ps_filter(Size == "S2") %>%
  comp_barplot(
    tax_level = "Order_Bin",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    facet_by = "Size",
    tax_transform_for_plot = "identity",
    sample_order=c( "TR08_B05_S2",  "TR11_B05_S2", "TR14_B05_S2"), 
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    #palette = custom_palette,
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) 


# ggsave("Abundance_S2_bins_20.pdf", plot = bin.abund.S2,
#        path = "plots/",
#        width = 10,
#        height = 5)

#--------------------

# S3,  2-0.2um filter fraction
bin.abund.S3= phylo_abund1  %>% ps_filter(Size == "S3") %>%
  comp_barplot(
    tax_level = "Order_Bin",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    facet_by = "Size",
    tax_transform_for_plot = "identity",
    sample_order=c( "TR08_B05_S3",  "TR11_B05_S3", "TR14_B05_S3"), 
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    #palette = custom_palette,
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) 


# ggsave("Abundance_S3_bins_20.pdf", plot = bin.abund.S3,
#        path = "plots/",
#        width = 10,
#        height = 5)
# 

#--------------------

# S4,  0.2-0.1um filter fraction
bin.abund.S4=  phylo_abund1  %>% ps_filter(Size == "S4") %>%
  comp_barplot(
    tax_level = "Order_Bin",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    facet_by = "Size",
    tax_transform_for_plot = "identity",
    sample_order=c("TR08_B05_S4",  "TR11_B05_S4", "TR12_B05_S4", "TR13_B05_S4","TR14_B05_S4"), 
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    #palette = custom_palette,
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) 


# ggsave("Abundance_S4_bins_20.pdf", plot = bin.abund.S4,
#        path = "plots/",
#        width = 10,
#        height = 5)


#--------------------

# S5,  0.1-0.025um filter fraction
bin.abund.S5=  phylo_abund1  %>% ps_filter(Size == "S5") %>%
  comp_barplot(
    tax_level = "Order_Bin",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    facet_by = "Size",
    tax_transform_for_plot = "identity",
    sample_order=c("TR08_B05_S5",  "TR11_B05_S5", "TR12_B05_S5", "TR13_B05_S5","TR14_B05_S5", "TR17_B05_S5"), 
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    #palette = custom_palette,
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) 


# ggsave("Abundance_S5_bins_20.pdf", plot = bin.abund.S5,
#        path = "plots/",
#        width = 10,
#        height = 5)
# 

