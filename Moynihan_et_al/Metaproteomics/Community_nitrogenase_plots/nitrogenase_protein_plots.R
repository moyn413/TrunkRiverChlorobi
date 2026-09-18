library(plyr)
library(dplyr)
library(microViz)
library(tidyverse)
library(phyloseq)
library(gridExtra)
library(RColorBrewer)


setwd("/path/")

N<-read.csv("nitrogenase_proteins.csv")

attach(N)

#Make phyloseq object for plotting

# Organize taxonomy and binID data
tax<-N %>% select(Gene:Tax_bin)
binID <- N %>% select(bin_gene)
tax<- cbind(binID, tax)

rownames(tax)<- tax$bin_gene
tax<- tax %>% select(-bin_gene)

#Protein abundance data
abund<-N%>%select(TR06_PB:TR17_M01_sup)
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



#=================================================================#
#==============Plots of nifH, nifD and nifK protein NSAF===================#
#=================================================================#
# Make subset objects that includes only bloom, mat and suspended biofilm ("web") from the B05 bloom sampling area
ps_BMW<-ps_gene%>%ps_filter(Type != "White") %>% ps_filter(Type2 != "Bloom_B01") %>% ps_filter(Type2 != "Supernatant_B01") %>% ps_filter(ID != "TR17_M01_DG")%>% ps_filter(ID != "TR17_M01_LG")

#------------------------------------------------------------------------
# Make plots of nifH, nifD, nifK--but first make color palette. 
# Need to remove genes to make color palette by order

tax_df <- as.data.frame(as.matrix(tax_table(ps_BMW)))  # Convert to dataframe for easier manipulation
tax_df <- tax_df[, !colnames(tax_df) %in% "Gene"]  # Remove Gene column (nifK, nifD,nifH)
tax_df <- tax_df[, !colnames(tax_df) %in% "Gene2"]  # Remove Gene column (nifK, nifD,nifH)
ps_nif_noG2<- ps_BMW
tax_table(ps_nif_noG2) <- tax_table(as.matrix(tax_df))

myPal2 <- tax_palette(data = ps_nif_noG2, rank = "Order", n = 25, pal = "brewerPlus")
tax_palette_plot(myPal2)


#------------------------------------------------------------------------
nifD <-ps_BMW %>% subset_taxa(Gene2 == "nifD")%>%ps_filter(Type2 == "Bloom_B05") %>% 
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    facet_by = "Type2",
    palette = myPal2,
    tax_transform_for_plot = "identity",
    # sample_order=c("TR06_PB", "TR08_B05", "TR09_B05","TR10_B05","TR11_B05", "TR12_B05","TR13_B05",
    #                "TR14_B05", "TR15_B05","TR16_B05","TR17_B05_MW","TR12_W_P","TR13_W_P","TR12_M01",
    #                "TR13_M01","TR14_M01", "TR16_M01", "TR17_M01", 
    #                "TR08_B05_sup","TR09_B05_sup","TR10_B05_sup","TR11_B05_sup", "TR13_B05_sup",
    #                "TR14_B05_sup","TR16_B05_sup","TR17_M01_sup"),
    sample_order=c("TR08_B05", "TR09_B05","TR10_B05","TR11_B05", "TR12_B05","TR13_B05","TR14_B05", "TR15_B05","TR16_B05","TR17_B05_MW"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  )  +  facet_wrap( facets = vars(Type2), scales = "free_x") + ggtitle("nifD")+ ylim(0,2)


#------------------------------------------------------------------------

nifK <-ps_BMW %>% subset_taxa(Gene2 == "nifK")%>%ps_filter(Type2 == "Bloom_B05") %>% 
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    facet_by = "Type2",
    palette = myPal2,
    tax_transform_for_plot = "identity",
    # sample_order=c("TR06_PB", "TR08_B05", "TR09_B05","TR10_B05","TR11_B05", "TR12_B05","TR13_B05",
    #                "TR14_B05", "TR15_B05","TR16_B05","TR17_B05_MW","TR12_W_P","TR13_W_P","TR12_M01",
    #                "TR13_M01","TR14_M01", "TR16_M01", "TR17_M01", 
    #                "TR08_B05_sup","TR09_B05_sup","TR10_B05_sup","TR11_B05_sup", "TR13_B05_sup",
    #                "TR14_B05_sup","TR16_B05_sup","TR17_M01_sup"),
    sample_order=c("TR08_B05", "TR09_B05","TR10_B05","TR11_B05", "TR12_B05","TR13_B05","TR14_B05", "TR15_B05","TR16_B05","TR17_B05_MW"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  )  +  facet_wrap( facets = vars(Type2), scales = "free_x") + ggtitle("nifK")+ ylim(0,2)

#------------------------------------------------------------------------
# Note that nifH sequences were clustered into I, II, or III based on phylogenies done in Geneious Prime against reference sequences. Only clusters I and III were present in this dataset. 
#------------------------------------------------------------------------
nifHI <-ps_BMW %>% subset_taxa(Gene2 == "nifH  Cluster I")%>%ps_filter(Type2 == "Bloom_B05") %>% 
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    facet_by = "Type2",
    palette = myPal2,
    tax_transform_for_plot = "identity",
    # sample_order=c("TR06_PB", "TR08_B05", "TR09_B05","TR10_B05","TR11_B05", "TR12_B05","TR13_B05",
    #                "TR14_B05", "TR15_B05","TR16_B05","TR17_B05_MW","TR12_W_P","TR13_W_P","TR12_M01",
    #                "TR13_M01","TR14_M01", "TR16_M01", "TR17_M01", 
    #                "TR08_B05_sup","TR09_B05_sup","TR10_B05_sup","TR11_B05_sup", "TR13_B05_sup",
    #                "TR14_B05_sup","TR16_B05_sup","TR17_M01_sup"),
    sample_order=c("TR08_B05", "TR09_B05","TR10_B05","TR11_B05", "TR12_B05","TR13_B05","TR14_B05", "TR15_B05","TR16_B05","TR17_B05_MW"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  )  +  facet_wrap( facets = vars(Type2), scales = "free_x") + ggtitle("nifH I")+ ylim(0,2)


#------------------------------------------------------------------------

nifHIII <-ps_BMW %>% subset_taxa(Gene2 == "nifH  Cluster III")%>%ps_filter(Type2 == "Bloom_B05") %>% 
  comp_barplot(
    tax_level = "Order",
    label = "SamplingDay", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    facet_by = "Type2",
    palette = myPal2,
    tax_transform_for_plot = "identity",
    # sample_order=c("TR06_PB", "TR08_B05", "TR09_B05","TR10_B05","TR11_B05", "TR12_B05","TR13_B05",
    #                "TR14_B05", "TR15_B05","TR16_B05","TR17_B05_MW","TR12_W_P","TR13_W_P","TR12_M01",
    #                "TR13_M01","TR14_M01", "TR16_M01", "TR17_M01", 
    #                "TR08_B05_sup","TR09_B05_sup","TR10_B05_sup","TR11_B05_sup", "TR13_B05_sup",
    #                "TR14_B05_sup","TR16_B05_sup","TR17_M01_sup"),
    sample_order=c("TR08_B05", "TR09_B05","TR10_B05","TR11_B05", "TR12_B05","TR13_B05","TR14_B05", "TR15_B05","TR16_B05","TR17_B05_MW"),
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  )  +  facet_wrap( facets = vars(Type2), scales = "free_x") + ggtitle("nifH III") + ylim(0,2)



x.order <-grid.arrange(nifHI, nifHIII, nifD, nifK , nrow = 1)

ggsave("nitrogenase_by_gene.pdf", plot = x.order, 
       path = "/path/", limitsize = FALSE,
       width = 30,
       height = 5)









