#=================================
# 16S rRNA processing of LoopSeq data using DADA2
#=================================

#=== This script follows the PacBio tutorial and modifications suggested on the DADA2 github issues page
#  https://benjjneb.github.io/LRASManuscript/LRASms_fecal.html
#  https://github.com/benjjneb/dada2/issues/991

#=== Clear all objects in workspace
rm(list = ls())

#=================================
#   Load Required Packages
#=================================

library(dada2)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(ShortRead)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(tibble)
library(readxl)
library(gridExtra)
library(readr)
library(stringr)
library(kableExtra)  

# =================================
# Switches & Variables 
# =================================

#=== If TRUE, testing takes only the first 4 fastq files
testing <- FALSE
# This turns on/off quality plot generation
do_plot_quality <- TRUE

#=== Primer sequences used by Loop for full length 16S rRNA sequencing
FWD <- "AGAGTTTGATCMTGGCTCAG"  #F27
REV <- "TACCTTGTTACGACTT"      #R1492

#=== Parameters for filterAndTrim
minLen=1000 
maxLen=2000
truncQ = 2        
maxEE=1
maxN = 0
minQ=3
OMEGA_A=1e-10

#=== parameters for removeBimeraDenovo
method_chimera = "consensus"


#=================================
#   File Name Parameters
#=================================

#=== Set working directory
setwd("/workingdirectorypath")

#=== Define path to dataset and to silva taxonomy database
dataset_path <- "/datasetpath"
database_file <- "~/silva_nr99_v138_wSpecies_train_set.fa.gz"


#=== Read contextual data file to build phyloseq. 
#=== Note that the sample names and order should match the names and order of samples in the DADA2 output table
sampleinfo <- read.csv("~/contextual_data.csv", header = TRUE)


#=== Define string to use as file identifier (modify as needed)
file_identifier = ".fq"
#=== Define string to use for name separation  (modify as needed)
file_name_separator = ".fq"     


#=================================
#   Set up file directories
#=================================

#=== Create functions to define file paths
path_dataset <- function(file_name) str_c(dataset_path, file_name)
file_dataset <- function(file_end)  str_c(dataset_path, "/dada2/", file_end)

#=== Define directory paths
fastq_dir <-    path_dataset("/fastq") # for raw sequences from illumina
qual_dir <-     path_dataset("/qual_pdf/")
dada2_dir <-    path_dataset("/dada2/")
blast_dir <-    path_dataset("/blast/")
cut_dir <- path_dataset("/cut")
filtered_dir <- path_dataset("/cut/filtered")

#=== Create directories
if(!dir.exists(qual_dir))dir.create(qual_dir) 
if(!dir.exists(dada2_dir))dir.create(dada2_dir) 
if(!dir.exists(blast_dir))dir.create(blast_dir)
if(!dir.exists(cut_dir))dir.create(cut_dir)
if(!dir.exists(fastq_dir))dir.create(fastq_dir)


# =================================
#    Get and organize file names
# =================================

#=== Get a list of all fastq files in the fastq directory and select those that have the "file_identifier"
fns <- sort(list.files(fastq_dir, full.names = TRUE)) 
if (testing) fns <- fns[1:4]
fns <- fns[str_detect(basename(fns), file_identifier)] 

#=== Extract sample names using the "file_name_separator"
sample.names <- str_split(basename(fns), pattern = file_name_separator, simplify = TRUE) 
sample.names <- sample.names[, 1] #select first column only

#=== Create an empty data frame and look at the number of sequences in each file
df <- data.frame()
# Loop through all the files 
for (i in 1:length(fns)) {
  # use the dada2 function fastq.geometry
  geom <- fastq.geometry(fns[i])
  # extract the information on number of sequences and file name
  df_one_row <- data.frame(n_seq = geom[1], file_name = basename(fns[i]))
  # add one line to data frame
  df <- bind_rows(df, df_one_row) }

#=== Save sequence information to .txt file
knitr::kable(df)
write.table(df, file = 'n_seq.txt', sep='\t', row.names = FALSE, na='',quote=FALSE)

# =================================
# Remove primers and put sequences in "cut" directory
# =================================

#=== Define paths for files to be put in cut directory
fns_cut <- file.path(cut_dir, basename(fns))

#=== Run removePrimers on all files. This will remove any sequences that do not have both the forward and reverse primers.
#    Only full length 16S rRNA sequences will have both forward and reverse primers. Note that for Loop sequencing, 
#    sequences that do not have primers can still be used, as long but non-full length (e.g. 1200bp) sequences may still be 
#    relevant and taxonomically informative 
#=== An alternative  would be to skip this step and keep all the sequences (with and without primers) and use only min and max length to filter sequences.
noprimers <- removePrimers(fns, fns_cut, primer.fwd=FWD, primer.rev=dada2:::rc(REV), max.mismatch = 2, orient=TRUE, verbose=TRUE)

#=== Look at length distribution
lens.fn <- lapply(fns_cut, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)

#===Extract sample names
get.sample.name <- function(fname) strsplit(basename(fname), file_name_separator)[[1]][1]
sample.names <- unname(sapply(fns_cut, get.sample.name))
head(sample.names)


# =================================
#   Plot quality after cutadapt
# =================================

if (do_plot_quality) {
  for(i in 1:length(fns_cut)) {
    print(str_c("i = ", i))
    p1 <- plotQualityProfile(fns_cut[i])
    # if (i <= 2) {print(p1)}
    p1_file <- paste0(qual_dir, unname(sapply(fns_cut[i], get.sample.name)),".qual.pdf")
    ggsave( plot=p1, filename= p1_file, device = "pdf", width = 15, height = 15, scale=1, units="cm")
    
    read_length <- data.frame(length = width(ShortRead::readFastq(fns_cut[i]))) # Read the fastQ file and get the length of all the reads...
    print(str_c("File before filtration", fns_cut[i], "- Read length min=", min(read_length$length),"max=", max(read_length$length), "mean=", mean(read_length$length, na.rm=TRUE),  sep=" "))
    
    print(paste("Finished with file", fns_cut[i], ". ", sep=""))
  }
}
  

# =================================
#   Filter and trim reads 
# =================================

#=== Filter and trim reads based on parameters define in beginning of script
fns_filt <- file.path(cut_dir, "filtered", basename(fns_cut))
out2 <- filterAndTrim(fns_cut, fns_filt, minQ=minQ, maxN = maxN, maxEE = maxEE, minLen=minLen, maxLen=maxLen,
                      truncQ = truncQ, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  
head(out2)

#=== In case not all files made it past this step, the file name list needs to be updated
fns_filt <- sort(list.files(filtered_dir, full.names = TRUE)) 
fns_filt <- fns_filt[str_detect(basename(fns_filt), file_identifier)] 


# =================================
#   Dereplicate, Error, and DADA2
# =================================

#=== Dereplicate reads 
derep <- derepFastq(fns_filt, verbose = FALSE)

#=== Name the derep-class objects by the sample names
names(derep) <- sample.names

#=== Learn and plot error rates
err <- learnErrors(derep, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=TRUE)
p <- plotErrors(err, nominalQ=TRUE)
p_file <- file_dataset("_LearnErrors.pdf")
ggsave(plot=p, filename= p_file, device = "pdf", 
        width = 15, height = 15, scale=1, units="cm")

#=== Run DADA2 on derep using err from above steps
dd2 <- dada(derep, err, BAND_SIZE=32, multithread=TRUE, OMEGA_A=OMEGA_A,DETECT_SINGLETONS=TRUE)
dd2[[1]]
saveRDS(dd2, file.path(dada2_dir, "TR2021_dd2.rds"))

# =================================
#    Track steps and write table 
# =================================   

track<-cbind(raw=noprimers[,1], primers=noprimers[,2], filtered=out2[,2], denoised=sapply(dd2, function(x) sum(x$denoised)))
write_tsv(data.frame(track), str_c(dada2_dir, "read_numbers_dada2.tsv"))

# =================================   
# Make sequence table 
# =================================   

st2 <- makeSequenceTable(dd2); dim(st2)

# =================================
#  Remove chimeras
# =================================  

seqtab.nochim <- removeBimeraDenovo(st2, method = method_chimera, multithread = TRUE, 
                                    verbose = TRUE)
# Compute % of non chimeras
paste0("% of non chimeras : ", sum(seqtab.nochim)/sum(st2) * 100)
paste0("total number of sequences : ", sum(seqtab.nochim))


# =================================
#  Assign taxonomy (this step takes a long time)
# =================================  
tax2 <- assignTaxonomy(seqs=seqtab.nochim, refFasta=database_file, multithread=TRUE, 
                       taxLevels = c("Domain","Phylum","Class","Order","Family","Genus", "Species"),
                       minBoot = 0, outputBootstraps = TRUE, verbose = TRUE)
saveRDS(tax2, str_c(dada2_dir, "TR2021_16S_Loop.taxa.rds"))

# =================================
#    Transforming and save the ASVs sequences and write output files
# =================================  
seqtab.nochim_trans <- as.data.frame(t(seqtab.nochim)) %>% rownames_to_column(var = "sequence") %>% 
  rowid_to_column(var = "OTUNumber") %>% mutate(OTUNumber = sprintf("otu%04d",OTUNumber)) %>% 
  mutate(sequence = str_replace_all(sequence, "(-|\\.)", ""))

df <- seqtab.nochim_trans
seq_out <- Biostrings::DNAStringSet(df$sequence)
names(seq_out) <- df$OTUNumber

#=== Write fasta file without taxonomy
Biostrings::writeXStringSet(seq_out, str_c(dada2_dir, "TR2021_16S_Loop_no_taxo.fasta"), 
                            compress = FALSE, width = 20000)


#=== Export the data 
taxa <- readRDS(str_c(dada2_dir, "TR2021_16S_Loop.taxa.rds"))
write.table(as_tibble(taxa), file = str_c(dada2_dir, "taxa.txt"), sep= "\t")
write_tsv(as_tibble(taxa$boot), file = str_c(dada2_dir, "taxa_boot.txt"))
write_tsv(as_tibble(seqtab.nochim), file = str_c(dada2_dir, "seqtab.txt"))


#=== Append taxonomy and boot to the sequence table 
taxa_tax <- as.data.frame(taxa$tax)
taxa_boot <- as.data.frame(taxa$boot) %>% rename_all(funs(str_c(., "_boot")))
seqtab.nochim_trans2 <- taxa_tax %>% bind_cols(taxa_boot) %>% bind_cols(seqtab.nochim_trans)


#=== Filter for 16S (bootstrap values can be filtered layer in phyloseq, so leaving 0 for now)
bootstrap_min <- 0
# Filter based on the bootstrap
seqtab.nochim_16S <- seqtab.nochim_trans2 %>% dplyr::filter(Domain_boot >= bootstrap_min)
# Create final table with bootstraps, sequences, abundances, and taxonomy. This file is very useful!
write_tsv(seqtab.nochim_16S, str_c(dada2_dir, "TR2021_16S_Loop_dada2.database.tsv"))

#=== Write file for BLAST
df <- seqtab.nochim_16S
seq_out <- Biostrings::DNAStringSet(df$sequence)
names(seq_out) <- str_c(df$OTUNumber, df$Domain, sep = "|")
Biostrings::writeXStringSet(seq_out, str_c(blast_dir, "TR2021_16S_ASV.fasta"), compress = FALSE, width = 20000)

# =================================
#   Reload sample names and files (if stopped above after taxonomic assignment)
# =================================

seqtab.nochim_16S <- read.table("~/TR2021_16S_Loop_dada2.database.tsv", header=TRUE, check.names = FALSE,sep = '\t')

#=== Extract sample names from seqtab
c1<-colnames(seqtab.nochim_16S)
c2<-c1[17:length(c1)]  #Modify if needed to select only sample names from header list
sample.names <- c2

# =================================
#    Make Phyloseq (see code above if need to reload and remake phyloseq)
# =================================

#Filter for 16S (remove bootstrap less than 90 in OTU file)
bootstrap_min <- 90
# Filter based on the bootstrap

seqtab.nochim_16S_boot <- seqtab.nochim_16S %>% dplyr::filter(Domain_boot >= bootstrap_min)

OTU <- seqtab.nochim_16S_boot %>%  remove_rownames(.) %>% column_to_rownames("OTUNumber") %>% select_if(is.numeric) %>% select(-contains("_boot")) %>% as.matrix() %>% otu_table(taxa_are_rows = TRUE)
TAX <- seqtab.nochim_16S_boot %>%  remove_rownames(.) %>% column_to_rownames("OTUNumber") %>% select(Domain:Species) %>% as.matrix() %>% tax_table()
seqs <- getSequences(seqtab.nochim_16S_boot$sequence)
names(seqs) <- seqtab.nochim_16S_boot$OTUNumber
sequences <- Biostrings::DNAStringSet(seqs)


#---------Organize contextual data---------------
samdf <- data.frame(sample_name = sample.names) 

# add additional contextual data categories as needed here
id <- sampleinfo$ID
name <-sampleinfo$Name
sampleday <-sampleinfo$SamplingDay
date <- sampleinfo$Date
type <- sampleinfo$SampleType
type2 <- sampleinfo$SampleType2
type3 <- sampleinfo$SampleType3


samdf <- data.frame(ID=id, Name=name, SamplingDay=sampleday, Date=date, SampleType=type, SampleType2=type2, SampleType3=type3)
rownames(samdf) <- sample.names


#---------save files---------------

ps_dada2 <- phyloseq(OTU, sample_data(samdf), TAX, sequences) #create phyloseq

saveRDS(ps_dada2, str_c(dada2_dir, "TR2021_Loop_16S_phyloseq.rds"))   #save phyloseq file

# Follow phyloseq tutorial for plotting (or see separate phyloseq script)
# https://joey711.github.io/phyloseq/import-data.html






