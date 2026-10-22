# Script used to process and analysis archaeal v4v5 amplicon data in the 2021 Trunk River GSB bloom dataset. 
# ------Load packages
library(plyr) 
library(dplyr)
library(ggplot2)
library(phyloseq)
library(decontam)
library(dada2)
library(kableExtra)
library(grid)
library(gridExtra)
library(tidyverse)
library(knitr)
library(microViz)

# -------------------------------------------#
# ------Load phyloseq object
# -------------------------------------------#
ps <- readRDS("~/TR2021_Arc_phyloseq.rds")   #read RDS phyloseq file


# -------------------------------------------#
# ------Address NAs in taxonomic ranking
# -------------------------------------------#

# Replace NA with "unspecified + [last identified taxa]"
# This prevents rows from being removed due to NAs in any downstream plotting
tax.clean <- data.frame(tax_table(ps))
for (i in 1:6){ tax.clean[,i] <- as.character(tax.clean[,i])}
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:6] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:6] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:6] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:6] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:6] <- family
  } 
}

tax_table(ps) <- as.matrix(tax.clean)


# -------------------------------------------#
# -----Remove contaminants and Chloroplasts
# -------------------------------------------#

# Using decontam tool, no contaminants found
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
contamdf.prev01 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.1)
table(contamdf.prev01$contaminant)
which(contamdf.prev01$contaminant)
#remove contaminants
ps.noncontam <- prune_taxa(!contamdf.prev01$contaminant, ps)

#Remove control samples
ps.samples <-  ps %>% subset_samples(Sample_or_Control != "Control") 

# Make separate phyloseqs for different sample types and document number of taxa.
# Note that sediment samples core are in this phyloseq object but not included in this study.
ps.all <- ps.samples
ps.bloom <-  ps.all %>% subset_samples(SampleType == "Bloom")

# Remove any ASVs that have an abundance of 0 
ps.all = prune_taxa(taxa_sums(ps.all)>=1, ps.all) #ntaxa =6478
ps.bloom = prune_taxa(taxa_sums(ps.bloom)>=1, ps.bloom) #ntaxa =2632

# Determine # of ASVs in mat, water, web for manuscript
ps.mat <- ps.all %>% subset_samples(SampleType == "Mat") 
ps.mat = prune_taxa(taxa_sums(ps.mat)>=1, ps.mat) #ntaxa =1382

ps.web <- ps.all %>% subset_samples(SampleType == "Web") 
ps.web = prune_taxa(taxa_sums(ps.web)>=1, ps.web) #ntaxa =971

ps.water <- ps.all %>% subset_samples(SampleType == "Water") 
ps.water = prune_taxa(taxa_sums(ps.water)>=1, ps.water) #ntaxa =2665


# -------------------------------------------#
## Inspect sample reads  
# -------------------------------------------#

# Show number of reads in each sample 
DNAreads <- ps.all %>% sample_sums(.) 
knitr::kable(DNAreads) %>% kable_styling(font_size = 8)
DNAreads.df <- as.data.frame(DNAreads)
write.table(DNAreads.df, file = "DNAreads.df.tsv", sep = "\t", row.names = FALSE)

# Remove samples with fewer than 10,000 reads 
# in the case of these samples, all samples have enough reads, but keeping this step as part of quality control
ps.all.LR <- prune_samples(sample_sums(ps.all)>=10000, ps.all)
ps.bloom.LR <- prune_samples(sample_sums(ps.bloom)>=10000, ps.bloom)


# -------------------------------------------#
## Make phyloseq with low abundant ASVs removed
# -------------------------------------------#
# Remove OTUs occurring fewer than 10x in entire dataset
ps.all.LO = prune_taxa(taxa_sums(ps.all.LR)>=10, ps.all.LR) #n=3951 (removes 2526 taxa)
ps.bloom.LO <- prune_taxa(taxa_sums(ps.bloom.LR)>=10, ps.bloom.LR) # n= 1660 (removes 972 taxa)

# -------------------------------------------#
# Noramalize
# -------------------------------------------#
ps_normalize_median <- function (ps, title) {
  ps_median = median(sample_sums(ps))
  cat(sprintf("\nThe median number of reads used for normalization of %s is  %.0f", title, ps_median))
  normalize_median = function(x, t=ps_median) (if(sum(x) > 0){ t * (x / sum(x))} else {x})
  ps = transform_sample_counts(ps, normalize_median)
  cat(str_c("\nPhyloseq ",title, "\n========== \n") )
  print(ps)
}

ps_abundant <- function(ps,contrib_min=0.01, title){
  total_per_sample <- max(sample_sums(ps))
  ps <- filter_taxa(ps, function(x) sum(x > total_per_sample*contrib_min) > 0, TRUE)
  ps <- ps_normalize_median(ps, title)
}


# Normalize phyloseq files with low abundant ASVs 
ps.all.LR.norm = ps_normalize_median(ps.all.LR, "Arc_all_LR_norm")
ps.bloom.LR.norm = ps_normalize_median(ps.bloom.LR, "Bloom_Arc_LR_norm")

# Normalize phyloseq file without low abundant ASVs
ps.all.norm = ps_normalize_median(ps.all.LO, "Arc_all_norm")
ps.bloom.norm = ps_normalize_median(ps.bloom.LO, "Bloom_LO_norm")

# Remove non-abundant taxa and normalize again
ps.all.abund = ps_abundant(ps.all.norm,contrib_min=0.01, "Arc_all_abund")  
ps.bloom.abund = ps_abundant(ps.bloom.norm,contrib_min=0.01, "Bloom_all_abund") 


#-------------------------------------------------------------------
# Bloom
#-------------------------------------------------------------------

bloom.bar = ps.bloom.norm %>% 
  comp_barplot(
    tax_level = "Genus",
    label = "Name", # name an alternative variable to label axis
    n_taxa = 16, # give more taxa unique colours
    # facet_by = "Size",
    # tax_transform_for_plot = "identity",
    sample_order=c("TR08_B05", "TR09_B05", "TR10_B05", "TR11_B05",
                   "TR12_B05","TR13_B05","TR14_B05","TR15_B05",
                   "TR16_B05","TR17_B05_MW"), 
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) 

# ggsave("TR_bloom_microviz_genus.pdf", plot = bloom.bar, 
#        path = "~/Plots/",
#        width = 8,
#        height = 5)



#-------------------------------------------------------------------
#Percent contribution calculations
#-------------------------------------------------------------------

#----------------------------------------------------------------
# Calculate % contribution of taxa of interest in main bloom period (excluding TR15 to TR17 timepoints)
earlybloom <- subset_samples(ps.bloom.LR.norm, Name!="TR15_B05" & Name!="TR16_B05" & Name!="TR17_B05")
methanomicro.bloom<-subset_taxa(earlybloom, Class == "Methanomicrobia")%>% taxa_sums(.)
methanosar.bloom<-subset_taxa(earlybloom, Class == "Methanosarcinia")%>% taxa_sums(.)
thermoplas.bloom<-subset_taxa(earlybloom, Class == "Thermoplasmata")%>% taxa_sums(.)
asgard.bloom<-subset_taxa(earlybloom, Phylum == "Asgardarchaeota")%>% taxa_sums(.)
aenigmarch.bloom<-subset_taxa(earlybloom, Phylum == "Aenigmarchaeota")%>% taxa_sums(.)

all <- earlybloom %>% taxa_sums(.)
x1 <-sum(methanomicro.bloom)/sum(all)*100;x1
y1 <-sum(methanosar.bloom)/sum(all)*100;y1
z1 <-sum(thermoplas.bloom)/sum(all)*100;z1
a1 <-sum(asgard.bloom)/sum(all)*100;a1
b1 <-sum(aenigmarch.bloom)/sum(all)*100;b1

# 15.83957 % Methanomicrobia
# 16.5315 % Methanosarcinia
# 31.30938 % Thermoplasmata
# 3.37602 % Asgardarchaeota
# 3.175031 % Aenigmarchaeota


#----------------------------------------------------------------
# Calculate % contribution of taxa of interest in later bloom (TR15 to TR17 timepoints))
latebloom <- subset_samples(ps.bloom.LR.norm, Name=="TR15_B05" | Name=="TR16_B05" | Name=="TR17_B05")
methanomicro.bloom<-subset_taxa(latebloom, Class == "Methanomicrobia")%>% taxa_sums(.)
methanosar.bloom<-subset_taxa(latebloom, Class == "Methanosarcinia")%>% taxa_sums(.)
thermoplas.bloom<-subset_taxa(latebloom, Class == "Thermoplasmata")%>% taxa_sums(.)
asgard.bloom<-subset_taxa(latebloom, Phylum == "Asgardarchaeota")%>% taxa_sums(.)
aenigmarch.bloom<-subset_taxa(latebloom, Phylum == "Aenigmarchaeota")%>% taxa_sums(.)


all <- latebloom %>% taxa_sums(.)
x1 <-sum(methanomicro.bloom)/sum(all)*100;x1
y1 <-sum(methanosar.bloom)/sum(all)*100;y1
z1 <-sum(thermoplas.bloom)/sum(all)*100;z1
a1 <-sum(asgard.bloom)/sum(all)*100;a1
b1 <-sum(aenigmarch.bloom)/sum(all)*100;b1

# 26.34602 % Methanomicrobia
# 4.496917 % Methanosarcinia
# 45.3318 % Thermoplasmata
# 5.24294 % Asgardarchaeota
# 1.324371 % Aenigmarchaeota


#----------------------------------------------------------------
# Calculate % Methogens and Thermoplasmata in mat
mat <- subset_samples(ps.all.LR.norm, SampleType=="Mat")
methanomicro.mat<-subset_taxa(mat, Class == "Methanomicrobia")%>% taxa_sums(.)
methanosar.mat<-subset_taxa(mat, Class == "Methanosarcinia")%>% taxa_sums(.)
thermoplas.mat<-subset_taxa(mat, Class == "Thermoplasmata")%>% taxa_sums(.)
all <- mat %>% taxa_sums(.)
x1 <-sum(methanomicro.mat)/sum(all)*100;x1
y1 <-sum(methanosar.mat)/sum(all)*100;y1
z1 <-sum(thermoplas.mat)/sum(all)*100;z1

# Mat
# 4.672193 % Methanomicrobia
# 8.514021 % Methanosarcinia
# 56.45112 % Thermoplasmata

#----------------------------------------------------------------
# Calculate % Methogens and Thermoplasmata in surface water (non bloom lagoon water)
water <- subset_samples(ps.all.LR.norm, SampleType=="Water")
methanomicro.water<-subset_taxa(water, Class == "Methanomicrobia")%>% taxa_sums(.)
methanosar.water<-subset_taxa(water, Class == "Methanosarcinia")%>% taxa_sums(.)
thermoplas.water<-subset_taxa(water, Class == "Thermoplasmata")%>% taxa_sums(.)
all <- mat %>% taxa_sums(.)
x1 <-sum(methanomicro.water)/sum(all)*100;x1
y1 <-sum(methanosar.water)/sum(all)*100;y1
z1 <-sum(thermoplas.water)/sum(all)*100;z1

# Water
# 10.16384 % Methanomicrobia
# 4.473737 % Methanosarcinia
# 44.26854 % Thermoplasmata


#----------------------------------------------------------------
# Calculate % Methogens and Thermoplasmata in web
web <- subset_samples(ps.all.LR.norm, SampleType=="Web")
methanomicro.web<-subset_taxa(web, Class == "Methanomicrobia")%>% taxa_sums(.)
methanosar.web<-subset_taxa(web, Class == "Methanosarcinia")%>% taxa_sums(.)
thermoplas.web<-subset_taxa(web, Class == "Thermoplasmata")%>% taxa_sums(.)
all <- mat %>% taxa_sums(.)
x1 <-sum(methanomicro.web)/sum(all)*100;x1
y1 <-sum(methanosar.web)/sum(all)*100;y1
z1 <-sum(thermoplas.web)/sum(all)*100;z1

# Water
# 1.286726 % Methanomicrobia
# 1.426385 % Methanosarcinia
# 30.77561 % Thermoplasmata



