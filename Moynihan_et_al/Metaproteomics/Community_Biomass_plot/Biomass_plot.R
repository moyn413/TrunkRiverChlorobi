library(plyr)
library(dplyr)
library(microViz)
library(tidyverse)
library(phyloseq)
library(RColorBrewer)

setwd("/Users/molly/Documents/GitHub/TrunkRiverChlorobi/Moynihan_et_al/Metaproteomics/Community_Biomass_plot")
biomass<-read.csv("Protein_Biomass.csv", header=TRUE)

rownames(biomass)<- biomass$Accession

# Organize taxonomy and binID data
tax<-biomass %>% select(1:7)%>%  mutate_if(is.character, as.factor)

#Make abundance table
abund <- biomass %>% select(9:45)

# Now I have my "otu" tables 
# And my taxonomy table, both with matching binIDs (aka otuIDs)
meta<-read.csv("Biomass_contextual_data.csv", header=TRUE)
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


# ggsave("Protein_biomass_B05_order_20_norm.pdf", plot = bloom ,
#        path = "Plots_Biomass/",
#        width = 10,
#        height = 5)



#--------------------
# Percent Contributions
#--------------------

# Subset data to Bloom05 bloom 
ps.B05 <- ps.norm %>% ps_filter(SampleType2 == "Bloom05")
# Merge taxa abundances based on taxonomic order
ps.B05.norm.Glommed = tax_glom(ps.B05, "Order")
# Remove any rows with 0s
ps.B05.norm.Glommed <- prune_taxa(taxa_sums(ps.B05.norm.Glommed) > 1, ps.B05.norm.Glommed)


#----------------------------------------------------------------
# Calculate % Chlorobiales in Bloom TR08
bloom <- subset_samples(ps.B05.norm.Glommed, Name=="TR08_B05")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
all <- bloom %>% taxa_sums(.)
x1 <-sum(chlorobi.bloom)/sum(all)*100;x1
# 71.82504 % Chloribales


#----------------------------------------------------------------
# Calculate % Chlorobiales in Bloom TR09
bloom <- subset_samples(ps.B05.norm.Glommed, Name=="TR09_B05")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
all <- bloom %>% taxa_sums(.)
x2 <-sum(chlorobi.bloom)/sum(all)*100;x2
# 77.95396 % Chlorobiales

#----------------------------------------------------------------
# Calculate % Chlorobiales in Bloom TR10
bloom <- subset_samples(ps.B05.norm.Glommed, Name=="TR10_B05")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
all <- bloom %>% taxa_sums(.)
x3 <-sum(chlorobi.bloom)/sum(all)*100;x3
# 81.17325 % Chlorobiales

#----------------------------------------------------------------
# Calculate % Chlorobiales in Bloom TR11
bloom <- subset_samples(ps.B05.norm.Glommed, Name=="TR11_B05")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
all <- bloom %>% taxa_sums(.)
x4 <-sum(chlorobi.bloom)/sum(all)*100;x4
# 85.42144 % Chlorobiales

#----------------------------------------------------------------
# Calculate % Chlorobiales in Bloom TR12
bloom <- subset_samples(ps.B05.norm.Glommed, Name=="TR12_B05")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
all <- bloom %>% taxa_sums(.)
x5 <-sum(chlorobi.bloom)/sum(all)*100;x5
# 90.27974 % Chlorobiales

#----------------------------------------------------------------
# Calculate % Chlorobiales in Bloom TR13
bloom <- subset_samples(ps.B05.norm.Glommed, Name=="TR13_B05")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
all <- bloom %>% taxa_sums(.)
x6 <-sum(chlorobi.bloom)/sum(all)*100;x6
# 90.38889 % Chlorobiales

#----------------------------------------------------------------
# Calculate % Chlorobiales in Bloom TR14
bloom <- subset_samples(ps.B05.norm.Glommed, Name=="TR14_B05")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
all <- bloom %>% taxa_sums(.)
x7 <-sum(chlorobi.bloom)/sum(all)*100;x7
# 76.8693 % Chlorobiales


#----------------------------------------------------------------
# Calculate % Chlorobiales in Bloom TR15
bloom <- subset_samples(ps.B05.norm.Glommed, Name=="TR15_B05")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
all <- bloom %>% taxa_sums(.)
x8 <-sum(chlorobi.bloom)/sum(all)*100;x8
# 60.72362 % Chlorobiales


#----------------------------------------------------------------
# Calculate % Chlorobiales in Bloom TR16
bloom <- subset_samples(ps.B05.norm.Glommed, Name=="TR16_B05")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
all <- bloom %>% taxa_sums(.)
x9 <-sum(chlorobi.bloom)/sum(all)*100;x9
# 3.986479 % Chlorobiales


#----------------------------------------------------------------
# Calculate % Chlorobiales in Bloom TR17
bloom <- subset_samples(ps.B05.norm.Glommed, Name=="TR17_B05")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
all <- bloom %>% taxa_sums(.)
x10 <-sum(chlorobi.bloom)/sum(all)*100;x10
# 2.578427 % Chlorobiales


#----------------------------------------------------------------
# Percent contributions of non-Chlorbiales
#----------------------------------------------------------------
# Calculate % Desulfovibrionales in Bloom B05
B05 <- subset_samples(ps.B05.norm.Glommed, SampleType2=="Bloom05")
desulf.B05<-subset_taxa(B05, Order == "Desulfovibrionales")%>% taxa_sums(.)
sums <- B05 %>% taxa_sums(.)
x11 <-sum(desulf.B05)/sum(sums)*100;x11
# 0.74 % Desulfovibrionales


# Calculate % Desulfobacterales in Bloom B05
B05 <- subset_samples(ps.B05.norm.Glommed, SampleType2=="Bloom05")
desulf.B05<-subset_taxa(B05, Order == "Desulfobacterales")%>% taxa_sums(.)
sums <- B05 %>% taxa_sums(.)
x12 <-sum(desulf.B05)/sum(sums)*100;x12
# 2.458712 % Desulfobacterales


# Calculate % Desulfurmonadales in Bloom B05
B05 <- subset_samples(ps.B05.norm.Glommed, SampleType2=="Bloom05")
desulf.B05<-subset_taxa(B05, Order == "Desulfuromonadales")%>% taxa_sums(.)
sums <- B05 %>% taxa_sums(.)
x13 <-sum(desulf.B05)/sum(sums)*100;x13
# 3.290117 % Desulfurmonadales


# Calculate % Cyanobacteriia in Bloom B05
B05 <- subset_samples(ps.B05.norm.Glommed, SampleType2=="Bloom05")
cyano.B05<-subset_taxa(B05, Class == "Cyanobacteriia")%>% taxa_sums(.)
sums <- B05 %>% taxa_sums(.)
x14 <-sum(cyano.B05)/sum(sums)*100;x14
# 8.124915 % Cyanobacteriia


# Calculate % Chromatiales in Bloom B05
B05 <- subset_samples(ps.B05.norm.Glommed, SampleType2=="Bloom05")
chrom.B05<-subset_taxa(B05, Order == "Chromatiales")%>% taxa_sums(.)
sums <- B05 %>% taxa_sums(.)
x15 <-sum(chrom.B05)/sum(sums)*100;x15
# 8.789057 % Chromatiales



# Calculate % Spirochaetales in Bloom B05
B05 <- subset_samples(ps.B05.norm.Glommed, SampleType2=="Bloom05")
spir.B05<-subset_taxa(B05, Order == "Spirochaetales")%>% taxa_sums(.)
sums <- B05 %>% taxa_sums(.)
x16 <-sum(spir.B05)/sum(sums)*100;x16
# 0.05180516 % Spirochaetales


#----------------------------------------------------------------
# Calculate % Semibin163 (GSB-TRL01) in Bloom B05
B05 <- subset_samples(ps.norm, SampleType2=="Bloom05")
BB163.bloom<-subset_taxa(B05, grepl("BB163", MergedBinID))%>% taxa_sums(.)
all <- B05 %>% taxa_sums(.)
y1 <-sum(BB163.bloom)/sum(all)*100;y1
# 59.54367 % Semibin163

# Calculate % metabat44 (GSB-TRL01)  in Bloom B05
B05 <- subset_samples(ps.norm, SampleType2=="Bloom05")
BA44.bloom<-subset_taxa(B05, grepl("BA44", MergedBinID))%>% taxa_sums(.)
all <- B05 %>% taxa_sums(.)
y2 <-sum(BA44.bloom)/sum(all)*100;y2
# 3.7 % metabat44

# Calculate % Semibin163 in peak bloom of Bloom B05
B05 <- subset_samples(ps.B05, SampleType2=="Bloom05")
B05.peak <- B05 %>% subset_samples(SamplingDay2!="TR15" & SamplingDay2!="TR16" & SamplingDay2!="TR17")
BB163.peakbloom<-subset_taxa(B05.peak, grepl("BB163", MergedBinID))%>% taxa_sums(.)
all <- B05.peak %>% taxa_sums(.)
y3 <-sum(BB163.peakbloom)/sum(all)*100;y3
# 76.2947 % Semibin163

# Calculate % Semibin163 in peak bloom of Bloom B05
B05 <- subset_samples(ps.B05, SampleType2=="Bloom05")
B05.peak <- B05 %>% subset_samples(SamplingDay2!="TR15" & SamplingDay2!="TR16" & SamplingDay2!="TR17")
BA44.peakbloom<-subset_taxa(B05.peak, grepl("BA44", MergedBinID))%>% taxa_sums(.)
all <- B05.peak %>% taxa_sums(.)
y3 <-sum(BA44.peakbloom)/sum(all)*100;y3
# 4.524791 % Semibin163




