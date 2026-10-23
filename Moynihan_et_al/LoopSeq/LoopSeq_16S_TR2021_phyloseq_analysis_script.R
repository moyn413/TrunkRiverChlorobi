# Script used to process and analysis LoopSeq full-length 16S rRNA sequencing data in the 2021 Trunk River GSB bloom dataset. 
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
# ------Load phyloseq object with LoopSeq data
# -------------------------------------------#
ps <- readRDS("/Users/molly/Documents/GitHub/TrunkRiverChlorobi/Moynihan_et_al/LoopSeq/TR2021_LoopSeq_16S_phyloseq.rds")   #read RDS phyloseq file

# -------------------------------------------#
# ------Address NAs in taxonomic ranking
# -------------------------------------------#


tax.clean <- data.frame(tax_table(ps))
for (i in 1:7){tax.clean[,i] <- as.character(tax.clean[,i])}
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    domain <- paste("Domain_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    genus <- paste("Genus_", tax.clean[i,6], sep = "")
    tax.clean[i, 7:7] <- genus
  } 
}

tax_table(ps) <- as.matrix(tax.clean)



# -------------------------------------------#
# -----Remove Chloroplasts & Mitochondria
# -------------------------------------------#
ps.noC <- subset_taxa(ps, Order != "Chloroplast"| is.na(Order)) # ntaxa=12515
ps.noC_noM <- subset_taxa(ps.noC, Family != "Mitochondria"| is.na(Family)) #ntaxa = 12463

#Make seperate chloroplast and mitochondira phyloseqs
ps.plastid <- subset_taxa(ps, Order == "Chloroplast")  #ntaxa = 317
ps.mito <- subset_taxa(ps, Family == "Mitochondria")   #ntaxa = 52

# -------------------------------------------#
# Determine # of ASVs in bloom, mat, water, suspended bio-film like matrix (aka "web") for manuscript
# -------------------------------------------#

# Make separate phyloseqs for different sample types and track number of taxa
ps.all <- ps.noC_noM
ps.all<- prune_taxa(taxa_sums(ps.all)>=1, ps.all) #ntaxa=12463

ps.bmw <-  ps.all %>% subset_samples(SampleType == "Bloom" | SampleType =="Web" | SampleType == "Mat")
ps.bmw<- prune_taxa(taxa_sums(ps.bmw)>=1, ps.bmw) #ntaxa= 5190

ps.bloom <-  ps.all %>% subset_samples(SampleType == "Bloom")
ps.bloom<- prune_taxa(taxa_sums(ps.bloom)>=1, ps.bloom) #ntaxa=3757

ps.water <-  ps.all%>% subset_samples(SampleType == "Water")
ps.water<- prune_taxa(taxa_sums(ps.water)>=1, ps.water) #ntaxa=2294

ps.web <-  ps.all%>% subset_samples(SampleType == "Web")
ps.web<- prune_taxa(taxa_sums(ps.web)>=1, ps.web) #ntaxa=747

ps.mat <-  ps.all%>% subset_samples(SampleType == "Mat")
ps.mat<- prune_taxa(taxa_sums(ps.mat)>=1, ps.mat) #ntaxa=1177


# -------------------------------------------#
## Inspect sample reads  (note that number of reads in LoopSeq is lower than in traditional short read amplicon sequencing)
# -------------------------------------------#

# Show number of reads in each sample 
DNAreads <- ps.all %>% sample_sums(.)
knitr::kable(DNAreads) %>% kable_styling(font_size = 8)

# Remove samples with fewer than 1000 reads (here, some sediment core samples not used in this paper but collected from the same system had low sample reads <1000)
ps.all.LR <- prune_samples(sample_sums(ps.all)>=1000, ps.all)
ps.bloom.LR <- prune_samples(sample_sums(ps.bloom)>=1000, ps.bloom)
ps.bmw.LR <- prune_samples(sample_sums(ps.bmw)>=1000, ps.bmw)


# -------------------------------------------#
## Make phyloseq with low abundant ASVs removed
# -------------------------------------------#
# Remove ASVs occurring fewer than 5x in entire dataset
ps.all.LO = prune_taxa(taxa_sums(ps.all.LR)>=5, ps.all.LR)
ps.bloom.LO <- prune_taxa(taxa_sums(ps.bloom.LR)>=5, ps.bloom.LR) 
ps.bmw.LO <- prune_taxa(taxa_sums(ps.bmw.LR)>=5, ps.bmw.LR)

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
ps.all.LR.norm = ps_normalize_median(ps.all.LR, "16S_all_LR_norm")
ps.bloom.LR.norm = ps_normalize_median(ps.bloom.LR, "Bloom_16S_LR_norm")
ps.bmw.LR.norm = ps_normalize_median(ps.bmw.LR, "BMW_16S_LR_norm")


# Normalize phyloseq file without low abundant ASVs
ps.all.norm = ps_normalize_median(ps.all.LO, "16S_all_norm")
ps.bloom.norm = ps_normalize_median(ps.bloom.LO, "Bloom_LO_norm")
ps.bmw.norm = ps_normalize_median(ps.bmw.LO, "Bloom_LO_norm")



# -------------------------------------------#
# MicroViz Abundance BarPlots
# -------------------------------------------#

custom_palette <- c("Chlorobiales" = "#FFE100",
                    "Campylobacterales"          = "#FFCC99" ,
                    "Chromatiales"          = "#C20088" ,
                    "Thiomicrospirales" = "#4C005C",
                    "Desulfuromonadaceae"  = "#5EF1F2" ,
                    "Rhodobacterales"          = "#a502f0" ,
                    "Cyanobacteriales"          = "#2BCE48" ,
                    "Bacteroidales" = "#00998F",
                    "Flavobacteriales"          = "#FFFF80" ,
                    "Sphingobacteriales"          = "#003380" ,
                    "Synechococcales"          = "#19A405" ,
                    "SAR11 clade"          = "#FF5000" ,
                    "Limnotrichales"= "#B2DF8A",
                    "Geopsychrobacteraceae" = "#fda613", 
                    "Spirochaetales"  = "#993F00",
                    "Vibrionales" = "#1F78B4", 
                    "Kiritimatiellales" = "#990000" ,
                    "JGI 0000069-P22"="#c2466b", 
                    "Bacteroidetes VC2.1 Bac22"= "#46c2a1", 
                    "Desulfobacterales"  = "#F0A3FF" ,
                    "Other genera" ="lightgrey" )

# Color palette above is set to match other amplicon plot. 
# Note that no prebloom timepoint available in LoopSeq dataset
 

ps.B05 <-ps.bmw.norm%>% subset_samples(SampleType3=="Bloom5"& SampleType=="Bloom")

bloom.bar= ps.B05 %>%
  comp_barplot(
    tax_level = "Order",
    label = "Name", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    sample_order=c("TR8-B05", "TR9-B05", "TR10-B05",
                   "TR11-B05",
                   "TR12-B05","TR13-B05","TR14-B05","TR15-B05",
                   "TR16-B05","TR17-B05-MW"), 
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    palette = custom_palette,
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)+
  ) 


 # ggsave("TR_Loop_bloom_microviz_20order_updated.pdf", plot = bloom.bar, 
 #        path = "~/Plots",
 #        width = 10,
 #        height = 5)





# -------------------------------------------#
# Percent Contribution
# -------------------------------------------#

#----------------------------------------------------------------
# Calculate % Chlorobi in Bloom B05--remove TR15, TR16, TR17 (late bloom)
bloom <- subset_samples(ps.bloom.LR.norm, SampleType3=="Bloom5")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(bloom), rownames(otu_table(bloom)) %in% c('otu0001'))
all <- bloom %>% taxa_sums(.)
x1 <-sum(chlorobi.bloom)/sum(all)*100;x1
y1 <-sum(asv1_subset)/sum(all)*100;y1
z1 <-sum(asv1_subset)/sum(chlorobi.bloom)*100;z1
# 69.13764  % Chloribales
# 53.2942 % otu0001 
# 77.0842 % otu0001 of Chlorobi



#----------------------------------------------------------------
# Calculate % Chlorobi in Bloom TR08
bloom <- subset_samples(ps.bloom.LR.norm, ID=="TR8-B05")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(bloom), rownames(otu_table(bloom)) %in% c('otu0001'))
all <- bloom %>% taxa_sums(.)
x2 <-sum(chlorobi.bloom)/sum(all)*100;x2
y2 <-sum(asv1_subset)/sum(all)*100;y2
z2 <-sum(asv1_subset)/sum(chlorobi.bloom)*100;z2

# 80.51999 % Chloribales
# 64.05357 % otu0001
# 79.5499 % otu0001 of Chlorobi



#----------------------------------------------------------------
# Calculate % Chlorobi in Bloom TR09
bloom <- subset_samples(ps.bloom.LR.norm, ID=="TR9-B05")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(bloom), rownames(otu_table(bloom)) %in% c('otu0001'))
all <- bloom %>% taxa_sums(.)
x3 <-sum(chlorobi.bloom)/sum(all)*100;x3
y3 <-sum(asv1_subset)/sum(all)*100;y3
z3 <-sum(asv1_subset)/sum(chlorobi.bloom)*100;z3

# 85.97938 % Chlorobiales
# 62.64292 % otu0001
# 72.85808 % otu0001 of Chlorobi

#----------------------------------------------------------------
# Calculate % Chlorobi in Bloom TR10
bloom <- subset_samples(ps.bloom.LR.norm, ID=="TR10-B05")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(bloom), rownames(otu_table(bloom)) %in% c('otu0001'))
all <- bloom %>% taxa_sums(.)
x4 <-sum(chlorobi.bloom)/sum(all)*100;x4
y4 <-sum(asv1_subset)/sum(all)*100;y4
z4 <-sum(asv1_subset)/sum(chlorobi.bloom)*100;z4

# 85.59307 % Chlorobiales
# 66.61472 % otu0001
# 77.82723 % otu0001 of Chlorobi


#----------------------------------------------------------------
# Calculate % Chlorobi in Bloom TR11
bloom <- subset_samples(ps.bloom.LR.norm, ID=="TR11-B05")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(bloom), rownames(otu_table(bloom)) %in% c('otu0001'))
all <- bloom %>% taxa_sums(.)
x5 <-sum(chlorobi.bloom)/sum(all)*100;x5
y5 <-sum(asv1_subset)/sum(all)*100;y5
z5 <-sum(asv1_subset)/sum(chlorobi.bloom)*100;z5
# 90.79507 % Chlorobiales
# 70.83987 % otu0001
# 78.02171 % otu0001 of Chlorobi



#----------------------------------------------------------------
# Calculate % Chlorobi in Bloom TR12
bloom <- subset_samples(ps.bloom.LR.norm, ID=="TR12-B05")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(bloom), rownames(otu_table(bloom)) %in% c('otu0001'))
all <- bloom %>% taxa_sums(.)
x6<-sum(chlorobi.bloom)/sum(all)*100;x6
y6 <-sum(asv1_subset)/sum(all)*100;y6
z6 <-sum(asv1_subset)/sum(chlorobi.bloom)*100;z6

# 96.56291 % Chlorobiales
# 75.07566 % otu0001
# 77.74793 % otu0001 of Chlorobi



#----------------------------------------------------------------
# Calculate % Chlorobi in Bloom TR13
bloom <- subset_samples(ps.bloom.LR.norm, ID=="TR13-B05")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(bloom), rownames(otu_table(bloom)) %in% c('otu0001'))
all <- bloom %>% taxa_sums(.)
x7 <-sum(chlorobi.bloom)/sum(all)*100;x7
y7 <-sum(asv1_subset)/sum(all)*100;y7
z7 <-sum(asv1_subset)/sum(chlorobi.bloom)*100;z7

# 92.09738 % Chlorobiales
# 71.68539 % otu0001
# 77.83652 % otu0001 of Chlorobi



#----------------------------------------------------------------
# Calculate % Chlorobi in Bloom TR14
bloom <- subset_samples(ps.bloom.LR.norm, Name=="TR14-B05")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(bloom), rownames(otu_table(bloom)) %in% c('otu0001'))
all <- bloom %>% taxa_sums(.)
x8 <-sum(chlorobi.bloom)/sum(all)*100;x8
y8 <-sum(asv1_subset)/sum(all)*100;y8
z8 <-sum(asv1_subset)/sum(chlorobi.bloom)*100;z8

# 91.39529 % Chlorobiales
# 71.00071 % otu0001
# 77.68531 % otu0001 of Chlorobi

# stopped here
#----------------------------------------------------------------
#----------------------------------------------------------------
#--------LATE BLOOM-------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
# Calculate % Chlorobi in Late Bloom TR15
latebloom <- subset_samples(ps.bloom.LR.norm, ID=="TR15-B05")
chlorobi.latebloom<-subset_taxa(latebloom, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(latebloom), rownames(otu_table(latebloom)) %in% c('otu0001'))
all <- latebloom %>% taxa_sums(.)
x9 <-sum(chlorobi.latebloom)/sum(all)*100;x9
y9 <-sum(asv1_subset)/sum(all)*100;y9
z9 <-sum(asv1_subset)/sum(chlorobi.latebloom)*100;z9

# 56.39269 % Chlorobiales
# 42.88976 % otu0001
# 76.05552 % otu0001 of Chlorobi

#----------------------------------------------------------------
# Calculate % Chlorobi in Late Bloom TR16
latebloom <- subset_samples(ps.bloom.LR.norm, ID=="TR16-B05")
chlorobi.latebloom<-subset_taxa(latebloom, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(latebloom), rownames(otu_table(latebloom)) %in% c('otu0001'))
all <- latebloom %>% taxa_sums(.)
x10 <-sum(chlorobi.latebloom)/sum(all)*100;x10
y10 <-sum(asv1_subset)/sum(all)*100;y10
z10 <-sum(asv1_subset)/sum(chlorobi.latebloom)*100;z10


# 9.605705 % Chlorobiales
# 6.040268 % otu0001
# 62.8821 % otu0001 of Chlorobi

#----------------------------------------------------------------
# Calculate % Chlorobi in Late Bloom TR17
latebloom <- subset_samples(ps.bloom.LR.norm, ID=="TR17-B05-MW")
chlorobi.latebloom<-subset_taxa(latebloom, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(latebloom), rownames(otu_table(latebloom)) %in% c('otu0001'))
all <- latebloom %>% taxa_sums(.)
x11 <-sum(chlorobi.latebloom)/sum(all)*100;x11
y11 <-sum(asv1_subset)/sum(all)*100;y11
z11 <-sum(asv1_subset)/sum(chlorobi.latebloom)*100;z11

# 2.434929 % Chlorobiales
# 2.099076 % otu0001
# 86.2069 % otu0001 of Chlorobi


# Non-Chlorobi percent contributions
#----------------------------------------------------------------
# Calculate % Chlorobi in Bloom B05-
bloom <- subset_samples(ps.bloom.LR.norm, SampleType3=="Bloom5")
taxa.bloom<-subset_taxa(bloom, Order == "Spirochaetales")%>% taxa_sums(.)
all <- bloom %>% taxa_sums(.)
spiro <-sum(taxa.bloom)/sum(all)*100;spiro
# 0.5921194  % Spirochaetales

#----------------------------------------------------------------
# Calculate % Thiomicrospirales in Bloom B05
bloom <- subset_samples(ps.bloom.LR.norm, SampleType3=="Bloom5")
taxa.bloom<-subset_taxa(bloom, Order == "Thiomicrospirales")%>% taxa_sums(.)
all <- bloom %>% taxa_sums(.)
thio <-sum(taxa.bloom)/sum(all)*100;thio
# 2.108702  % Thiomicrospirales 

#----------------------------------------------------------------
# Calculate % Chromatiales in Bloom B05
bloom <- subset_samples(ps.bloom.LR.norm, SampleType3=="Bloom5")
taxa.bloom<-subset_taxa(bloom, Order == "Chromatiales")%>% taxa_sums(.)
all <- bloom %>% taxa_sums(.)
psb <-sum(taxa.bloom)/sum(all)*100;psb
# 4.073675  % Chromatiales

#----------------------------------------------------------------
# Calculate % Desulfuromonadaceae in Bloom B05
bloom <- subset_samples(ps.bloom.LR.norm, SampleType3=="Bloom5")
taxa.bloom<-subset_taxa(bloom, Order == "Desulfuromonadaceae")%>% taxa_sums(.)
all <- bloom %>% taxa_sums(.)
desulf <-sum(taxa.bloom)/sum(all)*100;desulf
# 1.914328  % Desulfuromonadaceae 


#----------------------------------------------------------------
# Calculate % Bacteroidales in Bloom B05
bloom <- subset_samples(ps.bloom.LR.norm, SampleType3=="Bloom5")
taxa.bloom<-subset_taxa(bloom, Order == "Bacteroidales")%>% taxa_sums(.)
all <- bloom %>% taxa_sums(.)
bac <-sum(taxa.bloom)/sum(all)*100;bac
# 1.662898  % Bacteroidales 






