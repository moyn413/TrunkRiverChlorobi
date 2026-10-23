# Script used to process and analysis bacterial v4v5 amplicon data in the 2021 Trunk River GSB bloom dataset. 
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
# ------Load phyloseq with data from mutiple Trunk River sampling years
# -------------------------------------------#
#This phyloseq object contains data from multiple sampling years in Trunk River. 
# Only 2021 data were used in this study.
ps <- readRDS("~/TR2021-2023_16S_phyloseq.rds")   #read RDS phyloseq file

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
# ------Subset data for this study (2021 bloom  only)
# -------------------------------------------#

#Create object with only TR 2021 data + blanks (exclude 2021 sediment,  2022 data, and 2023 data from Trunk River; exlcude Namibia samples present on the same run)
ps.TR <-  ps %>% subset_samples(Experiment != "Namibia" & Experiment != "2021 Coring" & Experiment != "2022 Diel" & Experiment !="2022 Outflow" & Run != "2023")

# Samples TR11 and TR13 were reanalyzed on a separate run. 
# Reanalyzed samples included rerunning the same DNA samples and running re-extracted DNA samples. 
# After looking at bar plots of re-analyzed TR11 and TR13 samples, they look nearly identical. 
# I decided to keep the 'rerun' samples here in this study. 
ps.TR<- ps.TR %>% subset_samples(Name != "TR11_B05" & Name !="TR13_B05" & Name != "TR11_B05_Redo" & Name != "TR13_B05_Redo")  

# Save phyloseq object of 2021 data used in this study.
# This is the object that is shared on github. 
saveRDS(ps.TR, "~/Bv4v5_TR2021_phyloseq.rds")

# Read in file shared on github
ps.TR <- readRDS("~/Bv4v5_TR2021_phyloseq.rds")   #read RDS phyloseq file


# -------------------------------------------#
# -----Remove contaminants and Chloroplasts & Mitochondria
# -------------------------------------------#

# Run decontam in prevalence mode
sample_data(ps.TR)$is.neg <- sample_data(ps.TR)$samp_vs_control == "Control"
contamdf.prev01 <- isContaminant(ps.TR, method="prevalence", neg="is.neg", threshold=0.1)
table(contamdf.prev01$contaminant)
which(contamdf.prev01$contaminant) #8 contaminants identified
#remove contaminants
ps.TR.noncontam <- prune_taxa(!contamdf.prev01$contaminant, ps.TR)

#Remove control samples (n=20)
ps.TR.samples <-  ps.TR.noncontam %>% subset_samples(samp_vs_control != "Control")  #ntaxa=25665

# Remove any ASVs that have an abundance of 0 (the original phyloseq object had ASVs from other sampling years/projects that were on the same run, which were all 0 in this dataset since they were removed from the object)
ps.TR.samples = prune_taxa(taxa_sums(ps.TR.samples)>=1, ps.TR.samples) #ntaxa =3764

#Remove chloroplast sequences. (Subset taxa will remove NAs. In theory there should be no NAs here because of the tax.clean steps above, but it is good practice to keep this here so that data is not deleted.)
ps.noC <- subset_taxa(ps.TR.samples, Order != "Chloroplast"| is.na(Order)) # ntaxa=144
ps.noC_noM <- subset_taxa(ps.noC, Family != "Mitochondria"| is.na(Family)) #ntaxa = 54
ps.all <- ps.noC_noM #ntaxa=3566

#Make seperate chloroplast and mitochondira phyloseqs
ps.noncontam_plastid <- subset_taxa(ps.TR.samples, Order == "Chloroplast")  #ntaxa = 258
ps.noncontam_mito <- subset_taxa(ps.TR.samples, Family == "Mitochondria")   #ntaxa = 54

# -------------------------------------------#
# Determine # of ASVs in bloom, mat, water, web for manuscript
# -------------------------------------------#
ps.2021 <-  ps.all %>% subset_samples(Experiment == "2021 Timeseries") #ntaxa =3566

ps.2021.bloom <- ps.2021 %>% subset_samples(SampleType == "Bloom") 
ps.2021.bloom = prune_taxa(taxa_sums(ps.2021.bloom)>=1, ps.2021.bloom) #ntaxa =2207

ps.2021.mat <- ps.2021 %>% subset_samples(SampleType == "Mat") 
ps.2021.mat = prune_taxa(taxa_sums(ps.2021.mat)>=1, ps.2021.mat) #ntaxa =1365

ps.2021.web <- ps.2021 %>% subset_samples(SampleType == "Web") 
ps.2021.web = prune_taxa(taxa_sums(ps.2021.web)>=1, ps.2021.web) #ntaxa =767

ps.2021.water <- ps.2021 %>% subset_samples(SampleType == "Water") 
ps.2021.water = prune_taxa(taxa_sums(ps.2021.water)>=1, ps.2021.water) #ntaxa =1420


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
ps.2021.LR <- prune_samples(sample_sums(ps.2021)>=10000, ps.2021)

# -------------------------------------------#
## Make phyloseq with low abundant ASVs removed
# -------------------------------------------#
# Remove ASVs occurring fewer than 10x in entire dataset
ps.all.LO = prune_taxa(taxa_sums(ps.2021.LR)>=10, ps.2021.LR) #n=2363 (removes 1203 taxa)

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


# Normalize phyloseq files with all ASVs
ps.all.LR.norm = ps_normalize_median(ps.2021.LR, "16S_all_LR_norm")

# Normalize phyloseq file without low abundant ASVs
ps.all.norm = ps_normalize_median(ps.all.LO, "16S_all_norm")

# Normalize taxa that have a minimum contribution
# ps.all.abund = ps_abundant(ps.all.LR,contrib_min=0.005, "16S_all_abund")  


# -------------------------------------------#
# MicroViz Abundance BarPlots
# -------------------------------------------#

custom_palette <- c("Chlorobiales" = "#FFE100",
                    "Spirochaetales"  = "#993F00",
                    "Thiomicrospirales" = "#4C005C",
                    "Bacteroidales" = "#00998F",
                    "Enterobacterales"  = "#FF0010",
                    "SAR11 clade"          = "#FF5000" ,
                    "Chromatiales"          = "#C20088" ,
                    "Desulfuromonadales"          = "#5EF1F2" ,
                    "Flavobacteriales"          = "#FFFF80" ,
                    "Cyanobacteriales"          = "#2BCE48" ,
                    "Campylobacterales"          = "#FFCC99" ,
                    "Rhodobacterales"          = "#a502f0" ,
                    "Kiritimatiellales"          = "#990000" ,
                    "Desulfobacterales"          = "#F0A3FF" ,
                    "Sphingobacteriales"          = "#003380" ,
                    "Desulfovibrionales"          = "#426600" ,
                    "Fusobacteriales"          = "#E0FF66" ,
                    "Lachnospirales"          = "#94FFB5" ,
                    "Peptostreptococcales-Tissierellales"= "#808080" ,
                    "Synechococcales"          = "#19A405" ,
                    "Other genera" ="lightgrey" )



# Make object with just prebloom and bloom samples
ps.PB.B05 <-ps.all.norm%>% subset_samples(Name=="TR06_PB"| SampleType=="Bloom")

bloom.bar= ps.PB.B05 %>% 
  comp_barplot(
    tax_level = "Order",
    label = "Name", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    sample_order=c("TR06_PB", "TR08_B05", "TR09_B05", "TR10_B05",
                   "TR11_B05_rerun",
                   "TR12_B05","TR13_B05_rerun","TR14_B05","TR15_B05",
                   "TR16_B05","TR17_B05_MW"), 
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    palette = custom_palette,
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.9, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)+
  ) 


# ggsave("TR_Amplicon_prebloom_bloom_microviz_20order.pdf", plot = bloom.bar, 
#        path = "~/Plots/",
#        width = 10,
#        height = 5)




#-------------------
# Supplemental Bar Plots of Mat and Suspended Biofilm-like Matrix
#--------------------

mw.bar = ps.all.norm %>% ps_filter(SampleType == "Mat"| SampleType == "Web")%>%
  comp_barplot(
    tax_level = "Order",
    label = "Name", # name an alternative variable to label axis
    n_taxa = 10, # give more taxa unique colours
    facet_by = "SampleDetails",
    # tax_transform_for_plot = "identity",
    # sample_order=c("TR08_B05", "TR09_B05", "TR10_B05", "TR11_B05_rerun",
    #           "TR12_B05","TR13_B05_rerun","TR14_B05","TR15_B05",
    #          "TR16_B05","TR17_B05_MW"), 
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other genera", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) 

# ggsave("TR_mat_web_microviz_order_updated.pdf", plot = mw.bar, 
#        path = "~/Plots/",
#        width = 8,
#        height = 5)




#-------------------
#Percent contribution calculations
#-------------------

# Use normalized phyloseq object with all taxa (ps.all.LR.norm)
# Calculate percent Chlorobi, percent ASV1 and percent of Chlorobi that ASV1 represents

#----------------------------------------------------------------
# Calculate % Chlorobi in Bloom
bloom <- subset_samples(ps.all.LR.norm, SampleType2=="Bloom")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(bloom), rownames(otu_table(bloom)) %in% c('otu0001'))
all <- bloom %>% taxa_sums(.)
x1 <-sum(chlorobi.bloom)/sum(all)*100;x1
y1 <-sum(asv1_subset)/sum(all)*100;y1
z1 <-sum(asv1_subset)/sum(chlorobi.bloom)*100;z1
# 82.49687 % Chloribales
# 79.38266 % otu0001
# 96.22506 % of Chlorbiales are otu0001

#----------------------------------------------------------------
# Calculate % other abundant taxa in Bloom
bloom <- subset_samples(ps.all.LR.norm, SampleType2=="Bloom")
spiro.bloom<-subset_taxa(bloom, Order == "Spirochaetales")%>% taxa_sums(.)
psb.bloom<-subset_taxa(bloom, Order == "Chromatiales")%>% taxa_sums(.)
dsm.bloom<-subset_taxa(bloom, Order == "Desulfuromonadales")%>% taxa_sums(.)
bac.bloom<-subset_taxa(bloom, Class == "Bacteroidia")%>% taxa_sums(.)
cyano.bloom<-subset_taxa(bloom, Order == "Cyanobacteriales")%>% taxa_sums(.)

all <- bloom %>% taxa_sums(.)

spiro <-sum(spiro.bloom)/sum(all)*100;spiro
# 4.887691 % Spirochaetales
psb <-sum(psb.bloom)/sum(all)*100;psb
# 1.277024 % PSB
dsm <-sum(dsm.bloom)/sum(all)*100;dsm
# 2.240997 % Desulfuromonadales
bac <-sum(bac.bloom)/sum(all)*100;bac
# 2.766598 % Bacteroidia
cyano <-sum(cyano.bloom)/sum(all)*100;cyano
# 0.0622511 % Cyanobacteriales

#----------------------------------------------------------------
# Calculate % other abundant taxa in bloom that corresponds with web samples
bloom <- subset_samples(ps.all.LR.norm, Name=="TR12_B05" | Name == "TR13_B05_rerun")
spiro.bloom<-subset_taxa(bloom, Order == "Spirochaetales")%>% taxa_sums(.)
psb.bloom<-subset_taxa(bloom, Order == "Chromatiales")%>% taxa_sums(.)
dsm.bloom<-subset_taxa(bloom, Order == "Desulfuromonadales")%>% taxa_sums(.)
bac.bloom<-subset_taxa(bloom, Class == "Bacteroidia")%>% taxa_sums(.)
cyano.bloom<-subset_taxa(bloom, Order == "Cyanobacteriales")%>% taxa_sums(.)

all <- bloom %>% taxa_sums(.)

spiro <-sum(spiro.bloom)/sum(all)*100;spiro
# 2.249706 % Spirochaetales
psb <-sum(psb.bloom)/sum(all)*100;psb
# 0.3685184 % PSB
dsm <-sum(dsm.bloom)/sum(all)*100;dsm
# 2.644706 % Desulfuromonadales
bac <-sum(bac.bloom)/sum(all)*100;bac
# 0.793641 % Bacteroidia
cyano <-sum(cyano.bloom)/sum(all)*100;cyano
# 0.1191598 % Cyanobacteriales



#----------------------------------------------------------------
# Calculate % Chlorobi in Bloom TR08
bloom <- subset_samples(ps.all.LR.norm, Name=="TR08_B05")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(bloom), rownames(otu_table(bloom)) %in% c('otu0001'))
all <- bloom %>% taxa_sums(.)
x2 <-sum(chlorobi.bloom)/sum(all)*100;x2
y2 <-sum(asv1_subset)/sum(all)*100;y2
# 68.42582 % Chloribales
# 67.68236 % otu0001
z2 <-sum(asv1_subset)/sum(chlorobi.bloom)*100;z2
# 98.91348 of Chlorobi are otu0001

#----------------------------------------------------------------
# Calculate % Chlorobi in Bloom TR09
bloom <- subset_samples(ps.all.LR.norm, Name=="TR09_B05")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(bloom), rownames(otu_table(bloom)) %in% c('otu0001'))
all <- bloom %>% taxa_sums(.)
x3 <-sum(chlorobi.bloom)/sum(all)*100;x3
y3 <-sum(asv1_subset)/sum(all)*100;y3
# 87.66523 % Chlorobiales
# 87.2431 % otu0001
z3 <-sum(asv1_subset)/sum(chlorobi.bloom)*100;z3
# 99.51848 of Chlorobi are otu0001


#----------------------------------------------------------------
# Calculate % Chlorobi in Bloom TR10
bloom <- subset_samples(ps.all.LR.norm, Name=="TR10_B05")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(bloom), rownames(otu_table(bloom)) %in% c('otu0001'))
all <- bloom %>% taxa_sums(.)
x4 <-sum(chlorobi.bloom)/sum(all)*100;x4
y4 <-sum(asv1_subset)/sum(all)*100;y4
# 80.1559 % Chlorobiales
# 79.83556 % otu0001
z4 <-sum(asv1_subset)/sum(chlorobi.bloom)*100;z4
# 99.60035 of Chlorobi are otu0001



#----------------------------------------------------------------
# Calculate % Chlorobi in Bloom TR11
bloom <- subset_samples(ps.all.LR.norm, Name=="TR11_B05_rerun")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(bloom), rownames(otu_table(bloom)) %in% c('otu0001'))
all <- bloom %>% taxa_sums(.)
x5 <-sum(chlorobi.bloom)/sum(all)*100;x5
y5 <-sum(asv1_subset)/sum(all)*100;y5
# 84.13471 % Chlorobiales
# 75.12395 % otu0001
z5 <-sum(asv1_subset)/sum(chlorobi.bloom)*100;z5
# 89.29008 of Chlorobi are otu0001



#----------------------------------------------------------------
# Calculate % Chlorobi in Bloom TR12
bloom <- subset_samples(ps.all.LR.norm, Name=="TR12_B05")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(bloom), rownames(otu_table(bloom)) %in% c('otu0001'))
all <- bloom %>% taxa_sums(.)
x6<-sum(chlorobi.bloom)/sum(all)*100;x6
y6 <-sum(asv1_subset)/sum(all)*100;y6
# 96.07716 % Chlorobiales
# 95.42601 % otu0001
z6 <-sum(asv1_subset)/sum(chlorobi.bloom)*100;z6
# 99.60035 of Chlorobi are otu0001


#----------------------------------------------------------------
# Calculate % Chlorobi in Bloom TR13
bloom <- subset_samples(ps.all.LR.norm, Name=="TR13_B05_rerun")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(bloom), rownames(otu_table(bloom)) %in% c('otu0001'))
all <- bloom %>% taxa_sums(.)
x7 <-sum(chlorobi.bloom)/sum(all)*100;x7
y7 <-sum(asv1_subset)/sum(all)*100;y7
# 86.44797 % Chlorobiales
# 76.25938 % otu0001
z7 <-sum(asv1_subset)/sum(chlorobi.bloom)*100;z7
# 88.21419 of Chlorobi are otu0001


#----------------------------------------------------------------
# Calculate % Chlorobi in Bloom TR14
bloom <- subset_samples(ps.all.LR.norm, Name=="TR14_B05")
chlorobi.bloom<-subset_taxa(bloom, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(bloom), rownames(otu_table(bloom)) %in% c('otu0001'))
all <- bloom %>% taxa_sums(.)
x8 <-sum(chlorobi.bloom)/sum(all)*100;x8
y8 <-sum(asv1_subset)/sum(all)*100;y8
# 74.88532 % Chlorobiales
# 74.42035 % otu0001
z8 <-sum(asv1_subset)/sum(chlorobi.bloom)*100;z8
# 99.37909 of Chlorobi are otu0001


#----------------------------------------------------------------
# Calculate % Chlorobi in Late Bloom TR15
latebloom <- subset_samples(ps.all.LR.norm, Name=="TR15_B05")
chlorobi.latebloom<-subset_taxa(latebloom, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(latebloom), rownames(otu_table(latebloom)) %in% c('otu0001'))
all <- latebloom %>% taxa_sums(.)
x9 <-sum(chlorobi.latebloom)/sum(all)*100;x9
y9 <-sum(asv1_subset)/sum(all)*100;y9
# 37.83276 % Chlorobiales
# 37.20344 % otu0001

z9 <-sum(asv1_subset)/sum(chlorobi.latebloom)*100;z9
# 97.07265 of Chlorobi are otu0001



#----------------------------------------------------------------
# Calculate % Chlorobi in Late Bloom TR16
latebloom <- subset_samples(ps.all.LR.norm, Name=="TR16_B05")
chlorobi.latebloom<-subset_taxa(latebloom, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(latebloom), rownames(otu_table(latebloom)) %in% c('otu0001'))
all <- latebloom %>% taxa_sums(.)
x10 <-sum(chlorobi.latebloom)/sum(all)*100;x10
y10 <-sum(asv1_subset)/sum(all)*100;y10

# 20.94656 % Chlorobiales
# 18.4726 % otu0001

z10 <-sum(asv1_subset)/sum(chlorobi.latebloom)*100;z10
# 88.18919 of Chlorobi are otu0001

#----------------------------------------------------------------
# Calculate % Chlorobi in Late Bloom TR17
latebloom <- subset_samples(ps.all.LR.norm, Name=="TR17_B05")
chlorobi.latebloom<-subset_taxa(latebloom, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(latebloom), rownames(otu_table(latebloom)) %in% c('otu0001'))
all <- latebloom %>% taxa_sums(.)
x11 <-sum(chlorobi.latebloom)/sum(all)*100;x11
y11 <-sum(asv1_subset)/sum(all)*100;y11

# 3.918974 % Chlorobiales
# 3.368667 % otu0001

z11 <-sum(asv1_subset)/sum(chlorobi.latebloom)*100;z11
# 85.95787 of Chlorobi are otu0001


#----------------------------------------------------------------
#--------Mat-------------------------
#----------------------------------------------------------------
# Calculate % Chlorobi in all Early Mat
mat <- subset_samples(ps.all.LR.norm, SampleType2=="Mat")
chlorobi.mat<-subset_taxa(mat, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(mat), rownames(otu_table(mat)) %in% c('otu0001'))
all <- mat %>% taxa_sums(.)
x1 <-sum(chlorobi.mat)/sum(all)*100;x1
y1 <-sum(asv1_subset)/sum(all)*100;y1
# 87.59468 % Chloribales
# 84.69236 % otu0001

#----------------------------------------------------------------
# Calculate % other abundant taxa in Early mat
mat <- subset_samples(ps.all.LR.norm, SampleType2=="Mat")
spiro.mat<-subset_taxa(mat, Order == "Spirochaetales")%>% taxa_sums(.)
psb.mat<-subset_taxa(mat, Order == "Chromatiales")%>% taxa_sums(.)
dsm.mat<-subset_taxa(mat, Order == "Desulfuromonadales")%>% taxa_sums(.)
bac.mat<-subset_taxa(mat, Class == "Bacteroidia")%>% taxa_sums(.)
cyano.mat<-subset_taxa(mat, Order == "Cyanobacteriales")%>% taxa_sums(.)

all <- mat %>% taxa_sums(.)

spiro <-sum(spiro.mat)/sum(all)*100;spiro
# 1.038713 % Spirochaetales
psb <-sum(psb.mat)/sum(all)*100;psb
# 1.518079 % PSB
dsm <-sum(dsm.mat)/sum(all)*100;dsm
# 1.060297 % Desulfuromonadales
bac <-sum(bac.mat)/sum(all)*100;bac
# 1.183663 % Bacteroidia
cyano <-sum(cyano.mat)/sum(all)*100;cyano
# 5.567738 % Cyanobacteriales

#----------------------------------------------------------------
# Calculate % Chlorobi in all Late Mat
mat <- subset_samples(ps.all.LR.norm, SampleType2=="Mat-late")
chlorobi.mat<-subset_taxa(mat, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(mat), rownames(otu_table(mat)) %in% c('otu0001'))
all <- mat %>% taxa_sums(.)
x2 <-sum(chlorobi.mat)/sum(all)*100;x2
y2 <-sum(asv1_subset)/sum(all)*100;y2
# 6.319105 % Chloribales
# 3.536231 % otu0001

#----------------------------------------------------------------
# Calculate % other abundant taxa in Late mat
mat <- subset_samples(ps.all.LR.norm, SampleType2=="Mat-late")
spiro.mat<-subset_taxa(mat, Order == "Spirochaetales")%>% taxa_sums(.)
psb.mat<-subset_taxa(mat, Order == "Chromatiales")%>% taxa_sums(.)
dsm.mat<-subset_taxa(mat, Order == "Desulfuromonadales")%>% taxa_sums(.)
bac.mat<-subset_taxa(mat, Class == "Bacteroidia")%>% taxa_sums(.)
cyano.mat<-subset_taxa(mat, Order == "Cyanobacteriales")%>% taxa_sums(.)

all <- mat %>% taxa_sums(.)

spiro <-sum(spiro.mat)/sum(all)*100;spiro
# 1.510646 % Spirochaetales
psb <-sum(psb.mat)/sum(all)*100;psb
# 9.857747 % PSB
dsm <-sum(dsm.mat)/sum(all)*100;dsm
# 0.8403988 % Desulfuromonadales
bac <-sum(bac.mat)/sum(all)*100;bac
# 12.04203 % Bacteroidia
cyano <-sum(cyano.mat)/sum(all)*100;cyano
# 56.23872 % Cyanobacteriales


#----------------------------------------------------------------
#--------Suspended Biofilm-Like Matrix (aka "Web")---------------
#----------------------------------------------------------------
# Calculate % Chlorobi in all Web
web <- subset_samples(ps.all.LR.norm, SampleType2=="Web")
chlorobi.web<-subset_taxa(web, Order == "Chlorobiales")%>% taxa_sums(.)
asv1_subset <- subset(otu_table(web), rownames(otu_table(web)) %in% c('otu0001'))
all <- web %>% taxa_sums(.)
x1 <-sum(chlorobi.web)/sum(all)*100;x1
y1 <-sum(asv1_subset)/sum(all)*100;y1
# 85.12825 % Chloribales
# 82.56749 % otu0001

#----------------------------------------------------------------
# Calculate % other abundant taxa in Web
web <- subset_samples(ps.all.LR.norm, SampleType2=="Web")
spiro.web<-subset_taxa(web, Order == "Spirochaetales")%>% taxa_sums(.)
psb.web<-subset_taxa(web, Order == "Chromatiales")%>% taxa_sums(.)
dsm.web<-subset_taxa(web, Order == "Desulfuromonadales")%>% taxa_sums(.)
bac.web<-subset_taxa(web, Class == "Bacteroidia")%>% taxa_sums(.)
cyano.web<-subset_taxa(web, Order == "Cyanobacteriales")%>% taxa_sums(.)

all <- web %>% taxa_sums(.)

spiro <-sum(spiro.web)/sum(all)*100;spiro
# 0.4312601 % Spirochaetales
psb <-sum(psb.web)/sum(all)*100;psb
# 5.419174 % PSB
dsm <-sum(dsm.web)/sum(all)*100;dsm
# 0.9107486 % Desulfuromonadales
bac <-sum(bac.web)/sum(all)*100;bac
# 2.822383 % Bacteroidia
cyano <-sum(cyano.web)/sum(all)*100;cyano
# 0.1130753 % Cyanobacteriales


#----------------------------------------------------------------
# Calculate % other abundant taxa in bloom that correspond with web samples
bloom <- subset_samples(ps.all.LR.norm, SampleDetails=="Bloom5")
spiro.bloom<-subset_taxa(bloom, Order == "Spirochaetales")%>% taxa_sums(.)
psb.bloom<-subset_taxa(bloom, Order == "Chromatiales")%>% taxa_sums(.)
dsm.bloom<-subset_taxa(bloom, Order == "Desulfuromonadales")%>% taxa_sums(.)
dsulfb.bloom<-subset_taxa(bloom, Order == "Desulfobacterales")%>% taxa_sums(.)
dsulfbv.bloom<-subset_taxa(bloom, Order == "Desulfovibrionales")%>% taxa_sums(.)
cyano.bloom<-subset_taxa(bloom, Order == "Cyanobacteriales")%>% taxa_sums(.)

all <- bloom %>% taxa_sums(.)

spiro <-sum(spiro.bloom)/sum(all)*100;spiro
# 4.624718 % Spirochaetales
psb <-sum(psb.bloom)/sum(all)*100;psb
# 2.338124 % PSB
dsm <-sum(dsm.bloom)/sum(all)*100;dsm
# 2.043753 % Desulfuromonadales
desulfb <-sum(dsulfb.bloom)/sum(all)*100;desulfb
# 1.102853 % Desulfobacterales
dsulfbv <-sum(dsulfbv.bloom)/sum(all)*100;dsulfbv
# 0.8062449 % Desulfovibrionales
cyano <-sum(cyano.bloom)/sum(all)*100;cyano
# 1.929522 % Cyanobacteriales




