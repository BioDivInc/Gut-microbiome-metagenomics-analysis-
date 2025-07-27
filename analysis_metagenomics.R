# 0.0 Load dependencies --------

library(tidyverse) # for data manipulation, visualization, ...
library(ShortRead) # reading and examining of raw sequence reads
library(BiocManager) # provides a platform to install and manage packages for the e.g., analysis of genomic data
library(microbiome) # for microbiome composition visualization
library(dada2) # de-noise, trim, ..., pre-processing of raw sequencing data
library(rbiom) # for e.g., rarefaction performance
library(phyloseq) # further analysis of pre-processed sequencing data
library(Biostrings) # make sense of primer sequences
library(XML) # read, convert, ..., .xml meta data files
library(flextable) # to color-code data frames 
library(data.table) # enhanced data frame, useful in data wrangling
library(patchwork) # to align plots
library(DESeq2) # for differential abundance analysis
library(metagMisc) # to count unique taxa per group
library(ggvenn) # for venn diagrams
library(rstatix) # for statistical evaluations
library(vegan) # compute bray-curtis distance matrix and PERMANOVA
library(ggsignif) # add statistics results to plots
library(tibble) # to e.g., deframe() a df

# 1.0 Set path and import data -------------------------------------------------------------

# set path 
setwd("C:/Users/Felix/Desktop/R/MetaG/data")
path <- getwd()
list.files(path)

# use `sort` and `pattern` to ensure that forward and reverse reads are matched correctly
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))  # raw forward reads
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))  # raw reverse reads
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) # define raw sample names and get rid of information about forward or reverse reads


# 1.1 Cut primer sequences using `Biostrings`, `DECIPHER` and `cutadapt` --------

# forward primer (341F): CCTACGGG*N*GGC*W*GCAG
# reverse primer (805R): GACTAC*H**V*GGGTATCTAATCC
# N = any base, W = weak base (A or T), H = not G, V = not T
head(sread(readFastq(fnFs[1])))
head(sread(readFastq(fnFs[5])))
head(sread(readFastq(fnRs[1])))
head(sread(readFastq(fnRs[5])))

# set up primer removal using `Biostrings` to handle ambiguous bases and pass to `cutadapt` for cutting
# primer sequences used
Fs_primer <- Biostrings::DNAStringSet("CCTACGGGNGGCWGCAG")
Rs_primer <- Biostrings::DNAStringSet("GACTACHVGGGTATCTAATCC")

# disambiguate using `DECIPHER`
Fs_deciphered <- DECIPHER::Disambiguate(Fs_primer)
Rs_deciphered <- DECIPHER::Disambiguate(Rs_primer)

# get sequences out of the lists, convert to character and display
Fs_unlisted <- unlist(Fs_deciphered)
Rs_unlisted <- unlist(Rs_deciphered)
Fs_chr <- as.character(Fs_unlisted)
Rs_chr <- as.character(Rs_unlisted)
cat("Found", length(Fs_chr), "forward primer sequences:", Fs_chr)
cat("Found", length(Rs_chr), "reverse primer sequences:", Rs_chr)

# set path for `cutadapt`
cutadapt <- "cutadapt"  
system2(cutadapt, args = "--version")

# create new directory to store the clipped/cut files
cut_dir <- file.path(path, "cutadapt") # create new folder named "cutadapt"
if (!dir.exists(cut_dir)) dir.create(cut_dir) # if directory not yet exists, create it

# set up `cutadapt` and set arguments
fnFs_cut <- character(length(sample.names))
fnRs_cut <- character(length(sample.names))

# convert Fs_chr and Rs_chr to vector 
g <- as.vector(rbind("-g", paste0("^", Fs_chr)))
G <- as.vector(rbind("-G", paste0("^", Rs_chr)))

# loop through provided samples
for (i in seq_along(sample.names)){
  
  # naming of clipped files
  cutFs <- file.path(cut_dir, paste0(sample.names[i], "_F_cut.fastq.gz"))
  cutRs <- file.path(cut_dir, paste0(sample.names[i], "_R_cut.fastq.gz"))
  
  # define arguments
  cutadapt_args <- c(
    g, # forward anchored primer
    G, # reverse anchored primer
    "-o", cutFs, # output file for forward reads
    "-p", cutRs, # output file for reverse reads
    "--discard-untrimmed", # discard reads that do not contain expected primer
    "--error-rate=0.05", # decrease error rate to lower potential false matches
    "--cores=8", # use as many cores as your CPU has, to speed up computation
    fnFs[i], fnRs[i] # input files
  )
  
  # print out what's happening
  cat("Cutting sample", sample.names[i], "... \n")
  
  # run `cutadapt` with defined arguments
  system2(cutadapt, args = cutadapt_args)
}

# get cut files 
cut_fnFs <- sort(list.files(cut_dir, pattern = "_F_cut.fastq.gz", full.names = TRUE))  # raw forward reads
cut_fnRs <- sort(list.files(cut_dir, pattern = "_R_cut.fastq.gz", full.names = TRUE))  # raw reverse reads
names(cut_fnFs ) <- sample.names # prefix # define raw sample names and get rid of information about forward or reverse reads


# 1.2 Visualize quality profiles, `filterAndTrim()` of cut files  --------------

# visualize quality scores at each base position to check if trimming needs to be done
# discard reads with quality below quality score of set threshold
where_to_trim_Fs <- function(THRESHOLD, forward_reads){

  # define empty df
  df = data.frame()
  
  # loop through forward reads data and print out the respective annotated quality profile
  for (i in seq_along(forward_reads)){
    
    # create df with data from quality profile plots
    data_plot <- plotQualityProfile(forward_reads[i])$data
    
    # use mean, as the median is less affected by outliers
    avg_scores <- data_plot %>%
      group_by(Cycle) %>%
      summarise(avg_quality = sum(Score * Count) / sum(Count))
    
    # identify the cycles where the average quality (per cycle) < 30 and select first row
    if (any(avg_scores$avg_quality < THRESHOLD)){
      low_quality_cycles <- avg_scores %>%
        filter(avg_quality < THRESHOLD) %>%
        filter(row_number()==1)
      
      # store data in df
      updated_df <- data.frame(
        sample = sample.names[i], 
        Fs_cut = low_quality_cycles$Cycle
      )
      
      # bind together
      df <- rbind(df, updated_df)
      
      # print out what's happening
      cat("Suggested trimming of", sample.names[i], "based on set threshold of", THRESHOLD, "at", low_quality_cycles$Cycle, "cycles with a value of", round(low_quality_cycles$avg_quality, digits=2),".", "\n")
      
      # raw plot
      plot_raw <- plotQualityProfile(forward_reads[i])
      
      # annotate plot
      plot_annotated <- plot_raw +
        
        # insert red rectangle to display what section could be cut off
        geom_rect(data=low_quality_cycles, aes(xmin = Cycle, xmax = Inf, ymin = -Inf, ymax = Inf), inherit.aes = F, fill = "#ff4d3e", alpha=0.1)+
        
        # insert label with cycle information
        geom_label(
          data = low_quality_cycles,
          aes(x = Cycle, y = 10, label = Cycle),
          color = "white",
          fill = "#ff4d3e",
          label = low_quality_cycles$Cycle, size=4,
          show.legend = F)
      
    } else {
      
      low_quality_cycles <- avg_scores %>%
        slice_tail(n = 1)
      
        # store data in df
        updated_df <- data.frame(
          sample = sample.names[i], 
          Fs_cut = low_quality_cycles$Cycle
        )
        
        # bind together
        df <- rbind(df, updated_df)
        
        # print out what's happening
        cat("No trimming needed for", sample.names[i], "based on set threshold of", THRESHOLD,".\n")
        
        # raw plot
        plot_raw <- plotQualityProfile(forward_reads[i])
        
        plot_annotated <- plot_raw+
          geom_text(
            data = low_quality_cycles,
            aes(x = 60, y = 0),
            color = "#66c2a5",
            label = "Trimming:  not needed", size=3.9,
            show.legend = F)
    }
    print(plot_annotated)
  }
  # rename df
  trim_at_df_Fs <- df
  
  # check data frame for completeness
  return(trim_at_df_Fs)
}
where_to_trim_Rs <- function(THRESHOLD, reverse_reads){
  
  # define empty df
  df = data.frame()
  
  # loop through forward reads data and print out the respective annotated quality profile
  for (i in seq_along(reverse_reads)){
    
    # create df with data from quality profile plots
    data_plot <- plotQualityProfile(reverse_reads[i])$data
    
    # use mean, as the median is less affected by outliers
    avg_scores <- data_plot %>%
      group_by(Cycle) %>%
      summarise(avg_quality = sum(Score * Count) / sum(Count))
    
    # identify the cycles where the average quality (per cycle) < 30 and select first row
    if (any(avg_scores$avg_quality < THRESHOLD)){
      low_quality_cycles <- avg_scores %>%
        filter(avg_quality < THRESHOLD) %>%
        filter(row_number()==1)
      
      # store data in df
      updated_df <- data.frame(
        sample = sample.names[i], 
        Rs_cut = low_quality_cycles$Cycle
      )
      
      # bind together
      df <- rbind(df, updated_df)
      
      # print out what's happening
      cat("Suggested trimming of", sample.names[i], "based on set threshold of", THRESHOLD, "at", low_quality_cycles$Cycle, "cycles with a value of", round(low_quality_cycles$avg_quality, digits=2),".", "\n")
      
      # raw plot
      plot_raw <- plotQualityProfile(reverse_reads[i])
      
      # annotate plot
      plot_annotated <- plot_raw +
        
        # insert red rectangle to display what section could be cut off
        geom_rect(data=low_quality_cycles, aes(xmin = Cycle, xmax = Inf, ymin = -Inf, ymax = Inf), inherit.aes = F, fill = "#ff4d3e", alpha=0.1)+
        
        # insert label with cycle information
        geom_label(
          data = low_quality_cycles,
          aes(x = Cycle, y = 10, label = Cycle),
          color = "white",
          fill = "#ff4d3e",
          label = low_quality_cycles$Cycle, size=4,
          show.legend = F)
      
    } else {
      
      low_quality_cycles <- avg_scores %>%
        slice_tail(n = 1)
      
      # store data in df
      updated_df <- data.frame(
        sample = sample.names[i], 
        Rs_cut = low_quality_cycles$Cycle
      )
      
      # bind together
      df <- rbind(df, updated_df)
      
      # print out what's happening
      cat("No trimming needed for", sample.names[i], "based on set threshold of", THRESHOLD,".\n")
      
      # raw plot
      plot_raw <- plotQualityProfile(reverse_reads[i])
      
      plot_annotated <- plot_raw+
        geom_text(
          data = low_quality_cycles,
          aes(x = 60, y = 0),
          color = "#66c2a5",
          label = "Trimming:  not needed", size=3.9,
          show.legend = F)
    }
    print(plot_annotated)
  }
  # rename df
  trim_at_df_Rs <- df
  
  # check data frame for completeness
  return(trim_at_df_Rs)
}

trim_at_df_Fs <- where_to_trim_Fs(THRESHOLD = 27, forward_reads = cut_fnFs)
trim_at_df_Rs <- where_to_trim_Rs(THRESHOLD = 26, reverse_reads = cut_fnRs)

# merge the forward and reverse data frame
# define amplicon length
calc_amplicon_length <- function(forward_primer_binding_site, reverse_primer_binding_site, Fs_primer, Rs_primer){ # takes in Fs_primer and Rs_primer DNAStringSet objects from Biostrings
  
  # convert to numeric
  forward_primer_binding_site <- as.numeric(forward_primer_binding_site)
  reverse_primer_binding_site <- as.numeric(reverse_primer_binding_site)
  
  # if not already a Biostrig(DNAStringSet) object, convert
  Fs_primer <- Biostrings::DNAStringSet(Fs_primer)
  Rs_primer <- Biostrings::DNAStringSet(Rs_primer)
  
  # calculate the length of the amplicon with the trimmed of primers 
  amplicon_length <- (reverse_primer_binding_site-forward_primer_binding_site+1)-(Biostrings::width(Fs_primer)+Biostrings::width(Rs_primer))
  cat("The amplicon length is", amplicon_length, "bp.")
  return(amplicon_length)
}

amplicon_length <- calc_amplicon_length(forward_primer_binding_site="341", reverse_primer_binding_site = "805", Fs_primer="CCTACGGGNGGCWGCAG", Rs_primer="GACTACHVGGGTATCTAATCC")

# merge data frames with forward and reverse cut positions, add overlap
trim_at_df <- trim_at_df_Fs %>%
  right_join(trim_at_df_Rs) %>%
  mutate(overlap = Fs_cut+Rs_cut-amplicon_length)

# view data frame
trim_at_df

# color code data frame where overlap < 12 nt
trim_at_df_flex <- flextable(trim_at_df)
trim_at_df_flex <- color(trim_at_df_flex, i = ~ `overlap` < 12, j = "overlap", color="brown1", part = "body")
trim_at_df_flex

# place filtered files in filtered subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz")) # file path to filtered and trimmed forward reads
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz")) # file path to filtered and trimmed reverse reads
names(filtFs) <- sample.names # prefix

# empty list to store output from `filterAndTrim`
out_list <- vector("list", length = nrow(trim_at_df))

# loop through each provided sample and apply `filterAndTrim` function
for (i in 1:nrow(trim_at_df)) {
  truncLen_at <- c(trim_at_df$Fs_cut[i], trim_at_df$Rs_cut[i]) # takes one value from Fs and Rs at the same time --> c(Fs, Rs)
  cat("Trimming", sample.names[i], "at", truncLen_at, "cycles, respectively.", "\n") # print out sample name and where function will trim sequences
  
  out_list[[i]] <- filterAndTrim(cut_fnFs[i], filtFs[i], 
                                 cut_fnRs[i], filtRs[i],
                                 truncLen = truncLen_at,
                                 maxN = 0,
                                 maxEE = c(2, 2),
                                 truncQ = 2,
                                 rm.phix = TRUE,
                                 compress = TRUE,
                                 multithread = FALSE)
}

# summary of reads.in and reads.out after filtering and trimming sequences
out_df <- do.call(rbind, out_list) # rbinds each output from `filterAndTrim`
rownames(out_df) <- sample.names # adds sample names as row names
out_df <- as.data.frame(out_df) %>%
  dplyr::mutate("reads.out (%)" = (reads.out/reads.in)*100) # convert to percent
out_df

# sanity check if `filterAndTrim` worked as intended
SAMPLE = 5

# check forward trimming position via summary and quality profiling
sample_Fs <- readFastq(filtFs[SAMPLE])
plotQualityProfile(filtFs[SAMPLE])
trimmed_at_Fs <- Biostrings::width(sread(sample_Fs))
summary(trimmed_at_Fs)

# check reverse trimming position via summary and quality profiling
sample_Rs <- readFastq(filtRs[SAMPLE])
plotQualityProfile(filtRs[SAMPLE])
trimmed_at_Rs <- Biostrings::width(sread(sample_Rs))
summary(trimmed_at_Rs)


# 1.3 Learn errors, merge forward and reverse reads and remove chimeras --------------------------------------------------------------------

# learn error rates to compare observed and estimated errors; passed to `plotErrors` to visualize alignments
errF <- learnErrors(filtFs, multithread=TRUE) # learns errors from filtered forward files
errR <- learnErrors(filtRs, multithread=TRUE) # learns errors from filtered reverse files

# plot estimated and observed error rates vs. observed quality scores for each base position
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# decrease computation time by using dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# apply sample.names to dereplicated objects
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# takes error-corrected data and 
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# merge paired-end filtered (Fs, Rs) reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
mergers

# create sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # number of ASVs

# distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE) 
dim(seqtab.nochim) # number of ASVs left after getting rid of chimeras

sum(seqtab.nochim)/sum(seqtab) # 0-1; proportion of reads after chimera removal.


# 1.4 Pipeline summary ----------------------------------------------------

# track reads through the pre-processing pipeline
getN <- function(x) sum(getUniques(x))

# if single object, and NOT a list:
if (class(dadaFs)!="list"){
  track <- cbind(
    out_df,
    denoisedF = getN(dadaFs),
    denoisedR = getN(dadaRs),
    merged    = getN(mergers),
    nonchim   = rowSums(seqtab.nochim)
  )
  
  # print out object class
  cat("Class is a", capture.output(str(dadaFs))[1])
} else {
  # if list:
  track <- cbind(
    out_df,
    denoisedF = sapply(dadaFs, getN),
    denoisedR = sapply(dadaRs, getN),
    merged    = sapply(mergers, getN),
    nonchim   = rowSums(seqtab.nochim)
  )
  
  # print out object class
  cat("Class is a", capture.output(str(dadaFs))[1])
}

# define row and colnames
colnames(track) <- c("input", "filtered", "filtered (%)", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

# add merge rate to df and convert row names to column
track_df <- as.data.frame(track) %>%
  dplyr::mutate("merge_rate (%)" = (merged/filtered) * 100)

track_df <- setDT(track_df, keep.rownames = TRUE)
colnames(track_df)[1] <- "sample"
track_df

# print track and color-code merge rates with less than 75%
track_flex <- flextable(track_df)
track <- flextable::color(track_flex, i = ~ `merge_rate (%)`< 75, j = "merge_rate (%)", color = "brown1", part = "body")
track

# 1.5 Handle samples with low merge rate -----------------------------

# select samples with a merge rate of < 75%
low_mr_samples <- track_df %>%
  filter(`merge_rate (%)` < 75) %>%
  pull(sample) # select only sample column

# display
low_mr_samples

# if the length of low_mr_samples == 0 --> no samples found with mr < 75%:
if (length(low_mr_samples) == 0){
  seqtab.nochim_cleaned <- seqtab.nochim
} else { # if there are samples with mr < 75%:
  
  # update seqtab.nochim and delete rows corresponding to samples with mr < 75%
  seqtab.nochim_cleaned <- seqtab.nochim[!rownames(seqtab.nochim) %in% low_mr_samples,] # keep rows that are NOT in low_mr_samples 
  seqtab.nochim_cleaned
  
  # check `track` again with updated seqtab.nochim
  nonchim = rowSums(seqtab.nochim_cleaned)
  print(nonchim)
  length(nonchim)
}


# 1.6 Assign taxonomy -----------------------------------------------------

# assign taxonomy(Genus and/or Species)
taxa <- assignTaxonomy(seqtab.nochim_cleaned, "./silva_nr99_v138.2_toSpecies_trainset.fa.gz", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# 2.0 Import meta data --------------------------------------------------------

# select files with .xml extension
list_meta <- sort(list.files(pattern = ".xml"))

# create empty df to store metadata
meta_combined = data.frame()

# loop through data
for (i in list_meta){
  
  # import .xml files
  meta <- xmlParse(i)
  
  # convert to df
  df_meta <- xmlToDataFrame(nodes = getNodeSet(meta, "//VALUE"))
  
  # slice off uninteresting rows
  df_meta_adj <- df_meta %>%
    dplyr::slice(-c(1,2,3,5,9,10,11,12))
  
  # transpose and add column names
  t_df <- t(df_meta_adj)
  colnames(t_df ) <- c("Disease", "Sex", "Tissue", "Age")
  
  # append to combined df and thus update meta_combined
  meta_combined<- rbind(meta_combined, t_df)
  
}

# view combined df with updates row names
rownames(meta_combined) <- sample.names
meta_combined


# 3.0 Prepare data and data wrangling  -----------------------------------------------------

# combine multiple components in a single object
# otu_table(): OTU or ASV abundance data
# sample_data(): holds metadata
# tax_table(): stores taxonomy information
ps <- phyloseq(otu_table(seqtab.nochim_cleaned, taxa_are_rows=FALSE), 
               sample_data(meta_combined), 
               tax_table(taxa))
ps

# save file if it does not exist, else load it
if (!file.exists("./ps.rds")){
  print("File does not yet exist. Saving it...")
  saveRDS(ps, "ps.rds")
} else {
  ps <- readRDS("./ps.rds")
  print("File loaded successfully.")
}

# store in refseq slot in phyloseq object, rename sequences to ASVs
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# select top 20 abundant taxa for relative abundance plots
# top 20 abundant taxa, unnormalized
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- prune_taxa(top20, ps)

# subset Control
ps.top20_C <- ps.top20 %>%
  microbiome::transform(transform = "compositional") %>%
  aggregate_taxa(level = "Genus") %>%
  subset_samples(Disease == "Control")
ps.top20_C

# subset Parkinson
ps.top20_PD <- ps.top20 %>%
  microbiome::transform(transform = "compositional") %>%
  aggregate_taxa(level = "Genus") %>%
  subset_samples(Disease == "Parkinson disease")
ps.top20_PD

# 3.1 Visual analysis: Alpha diversity -----------------------------------------------------

# plot alpha-diversity and compare groups
alpha <- plot_richness(ps, x="Disease", measures=c("Chao1", "Shannon", "InvSimpson"))
alpha

# check layers and exclude geom_point() and geom_errorbar() to clean up the appearance
alpha$layers
alpha$layers <- alpha$layers[c(-1, -2)]

# change index names by indexing into data -> variables
levels(alpha$data[[6]]) <- c("Chao-1 Index", "Shannon Index", "Inverse Simpson Index")
levels(alpha$data[[6]])

# create df with annotation information
annotation_df <- data.frame(
  variable = levels(alpha$data[[6]]), # column name of indexes in alpha$data
  start = 1, # can also type in "Control"
  end = 2, # can also type in "Parkinson disease"
  y = c(800, 5, 60), # y-position
  label = c("NS", "NS", "NS") # label
)

# display df
annotation_df

# plot alpha diversity with added boxplots
alpha_diversity <- alpha + 
  geom_violin(aes(fill=Disease), alpha = 0.8) + 
  scale_fill_manual(values = c("cornflowerblue", "darkorange"))+
  geom_boxplot(width=0.1, linewidth = 0.6, colour = "black")+
  geom_point(position = position_dodge(width = 0.75)) + 
  theme_classic()+
  ylab("Alpha diversity measure")+
  xlab("Group")+
  geom_signif(
    data = annotation_df,
    aes(xmin = start, xmax = end, annotations = label, y_position = y),
    textsize = 5, size = 1, tip_length = 0, 
    manual = T
  ) +
  theme(
    axis.text.x = element_text(colour = "black", size = 13, face = "plain"),
    axis.text.y = element_text(colour = "black", size = 13, face = "plain"),
    axis.title.x = element_text(colour = "black", size = 16, face="bold"),
    axis.title.y = element_text(colour = "black", size = 16, face = "bold"),
    strip.text.x = element_text(colour = "black", size = 14, face = "bold"),
    legend.title = element_text(colour = "black", size = 14, face = "plain"),
    legend.text = element_text(colour = "black", size = 10, face = "plain"),
    line = element_line(linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    legend.position = "none"
  )

alpha_diversity

# save plot
ggsave("alpha_diversity.png", units="cm", width=35, height=25, dpi=800)
dev.off()


# 3.2 Visual analysis: Unique genera per sample and group -----------------

# get genera 
counts <- phyloseq_ntaxa_by_tax(ps, TaxRank = "Genus")
counts

# plot unique genera per and across sample(s)
counts_unique_sample <- counts %>%
  group_by(Sample, Disease) %>% # group by sample
  summarise(unique_genera = n_distinct(Genus)) %>% # unique number of genera per sample
  ungroup()

# display df
counts_unique_sample 

samples_control <- ggplot(subset(counts_unique_sample, Disease == "Control"), aes(x=Sample, y=unique_genera, fill = Disease))+
  geom_hline(yintercept = subset(counts_unique_sample, Disease == "Control") %>% summarise(mean(unique_genera)) %>% pull(), linewidth = 1, linetype = "dashed", color = "black", alpha = 0.8)+
  geom_bar(stat = "identity", color = "black", linewidth = 1, alpha = 0.8)+
  facet_grid(~Disease, scales="free")+
  scale_fill_manual(values = c("cornflowerblue"))+
  ylab("Number of unique genera")+
  theme_classic()+
  labs(tag = "A")+
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black", size = 13, face = "plain", angle = 90, vjust = .4),
        axis.text.y = element_text(colour = "black", size = 13, face = "plain"),
        axis.title.x = element_text(colour = "black", size = 16, face = "bold"),
        axis.title.y = element_text(colour = "black", size = 16, face = "bold"),
        strip.text.x = element_blank(),
        line = element_line(linewidth = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        plot.tag = element_text(hjust=0, size=25, color = "black"))

samples_control 

samples_PD <- ggplot(subset(counts_unique_sample, Disease == "Parkinson disease"), aes(x=Sample, y=unique_genera, fill = Disease))+
  geom_hline(yintercept = subset(counts_unique_sample, Disease == "Parkinson disease") %>% summarise(mean(unique_genera)) %>% pull(), linewidth = 1, linetype = "dashed", color = "black", alpha = 0.8)+
  geom_bar(stat = "identity", color = "black", linewidth = 1, alpha = 0.8)+
  facet_grid(~Disease, scales="free")+
  scale_fill_manual(values = "darkorange1")+
  theme_classic()+
  ylab("Number of unique genera")+
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black", size = 13, face = "plain", angle = 90, vjust = .4),
        axis.text.y = element_text(colour = "black", size = 13, face = "plain"),
        axis.title.x = element_text(colour = "black", size = 16, face = "bold"),
        axis.title.y = element_text(colour = "black", size = 16, face = "bold"),
        strip.text.x = element_blank(),
        line = element_line(linewidth = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

samples_PD

unique_genera_sum <- ggplot(counts_unique_sample, aes(x=Disease, y=unique_genera))+
  geom_violin(aes(fill=Disease), linewidth = 1, colour = "black", alpha = 0.8) + 
  scale_fill_manual(values = c("cornflowerblue", "darkorange"))+
  geom_boxplot(width=0.1, linewidth = 0.6, colour = "black")+
  geom_point(position = position_dodge(width = 0.75)) + 
  theme_classic()+
  ylab("Number of unique genera")+
  labs(fill="Group")+ # new legend title
  xlab("Group")+
  labs(tag = "B")+
  geom_signif(
    y_position = 160, xmin = 1, xmax = 2,
    annotation = c("NS"), tip_length = 0,
    size = 1,
    textsize = 5
  ) +
  theme(
    axis.text.x = element_text(colour = "black", size = 13, face = "plain", vjust = -15),
    axis.text.y = element_text(colour = "black", size = 13, face = "plain"),
    axis.title.x = element_text(colour = "black", size = 16, face="bold"),
    axis.title.y = element_text(colour = "black", size = 16, face = "bold"),
    strip.text.x = element_text(colour = "black", size = 14, face = "bold"),
    legend.title = element_text(colour = "black", size = 14, face = "plain"),
    legend.text = element_text(colour = "black", size = 10, face = "plain"),
    line = element_line(linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    plot.tag = element_text(hjust=0, size=25, color = "black")
  )

unique_genera_sum

# combine plots
per_sample <- (samples_control+samples_PD+unique_genera_sum)+plot_layout(guides = "collect", axes = "collect")
per_sample

# plot absolute number of unique taxa per group
# count distinct genera across samples per group
counts_unique_group <- counts %>%
  group_by(Disease) %>%
  summarise(unique_genera = n_distinct(Genus))

# display df
counts_unique_group

# plot difference using geom_bar()
unique_taxa_per_group <- ggplot(counts_unique_group , aes(x=Disease, y=unique_genera, fill = Disease))+
  geom_bar(stat = "identity", width = 0.8, colour = "black", linewidth = 1, alpha = 0.8)+
  scale_fill_manual(values = c("cornflowerblue", "darkorange"))+
  theme_classic()+
  ylab("Number of unique genera")+
  xlab("Group")+
  labs(tag = "C")+
  theme(
    axis.text.x = element_text(colour = "black", size = 13, face = "plain"),
    axis.text.y = element_text(colour = "black", size = 13, face = "plain"),
    axis.title.x = element_text(colour = "black", size = 16, face="bold"),
    axis.title.y = element_text(colour = "black", size = 16, face = "bold"),
    strip.text.x = element_text(colour = "black", size = 14, face = "bold"),
    legend.title = element_text(colour = "black", size = 14, face = "plain"),
    legend.text = element_text(colour = "black", size = 10, face = "plain"),
    line = element_line(linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    legend.position="none",
    plot.tag = element_text(hjust=0, size=25, color = "black")
  )

unique_taxa_per_group

# plot Venn diagram to explore putative unique genera depending on health 
# convert phyloseq object to data frame, filter out genera with an abundance of 0 as well as NAs
ps_df <- psmelt(ps) %>%
  filter(Abundance > 0) %>%
  filter(!is.na(Genus))

# select distinct 
ps_df_selected <- ps_df %>%
  distinct(Disease, Genus) %>%      
  group_by(Disease) %>%
  summarise(genera = list(unique(Genus))) %>%
  deframe()

# display df
ps_df_selected

# convert filtered data frame to list
control_genera <- ps_df_selected$`Control`
parkinson_genera <- ps_df_selected$`Parkinson disease`

# unique genera found in control
unique_genera_control <- setdiff(control_genera, parkinson_genera)
unique_genera_control

# unique genera found in parkinson disease samples
unique_genera_pd <- setdiff(parkinson_genera, control_genera)
unique_genera_pd

# plot venn diagram
venn_plot <- ggvenn::ggvenn(list(
  Control = control_genera,
  Parkinson = parkinson_genera),
  digits = 0,
  auto_scale = TRUE,
  fill_color = c("cornflowerblue", "darkorange1"),
  fill_alpha = 0.8,
  stroke_size = 1,
  set_name_size = 0,
  text_color = "black",
  text_size = 4.5)+
  labs(tag = "D")+
  theme_classic()+
  ylab("Difference in unique genera")+
  xlab("Group")+
  xlim(c(-4,4))+ # expand axis limits
  scale_x_continuous(breaks = c(-0.6, 0.6), labels = c("Control", "Parkinson disease")) + # add x-axis text
  theme(
    plot.tag = element_text(size=25, color = "black"),
    axis.title.y = element_text(colour = "black", size = 16, angle = 90, vjust = 2.5, face="bold"),
    axis.title.x = element_text(colour = "black", size = 16, face="bold"),
    axis.text.x = element_text(colour = "black", size = 13, face = "plain"),
    axis.text.y = element_blank(),
    axis.line = element_line(linewidth=1), # add axes lines to plot
    axis.ticks.x = element_line(linewidth = 1, color = "black"), # add ticks
    axis.ticks.y = element_blank(),
    axis.ticks.length = unit(0.1, "cm") # adjust tick length
  )

venn_plot

# combine plots
per_group <- (unique_taxa_per_group+venn_plot) + plot_layout(axes = "collect", guides = "collect")
per_group

# 3.3 Visual analysis: Microbiome composition -----------------------------

# plot composition per sample for Control
composition_per_sample_C <- plot_composition(ps.top20_C, 
                                           otu.sort = "abundance",
                                           group_by = "Disease")+
  geom_col(width = 1)+
  scale_fill_viridis_d("Genus")+
  theme_classic()+
  labs(tag = "A")+
  theme(
    axis.text.x = element_text(colour = "black", size = 13, face = "plain", angle = 90, vjust = .4),
    axis.text.y = element_text(colour = "black", size = 13, face = "plain"),
    axis.title.x = element_text(colour = "black", size = 16, face="bold"),
    axis.title.y = element_text(colour = "black", size = 16, face = "bold"),
    strip.text.x = element_text(colour = "black", size = 14, face = "bold"),
    legend.title = element_text(colour = "black", size = 14, face = "plain"),
    legend.text = element_text(colour = "black", size = 10, face = "plain"),
    line = element_line(linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position = "none",
    plot.tag = element_text(hjust=0, size=25, color = "black"),
    axis.ticks = element_line(colour = "black")
  )
composition_per_sample_C


# plot composition per group for Control
composition_per_group_C <- plot_composition(ps.top20_C, 
                                otu.sort = "abundance",
                                average_by = "Disease")+
  geom_col(width = 0.8)+
  scale_fill_viridis_d("Genus")+
  theme_classic()+
  labs(tag = "B")+
  scale_x_discrete(labels= "Average abundance")+
  theme(
    axis.text.x = element_text(colour = "black", size = 13, face = "plain", vjust = -15),
    axis.text.y = element_text(colour = "black", size = 13, face = "plain"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(colour = "black", size = 16, face = "bold"),
    strip.text.x = element_text(colour = "black", size = 14, face = "bold"),
    legend.title = element_text(colour = "black", size = 14, face = "plain"),
    legend.text = element_text(colour = "black", size = 10, face = "plain"),
    line = element_line(linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    plot.tag = element_text(hjust=0, size=25, color = "black"),
    axis.ticks = element_line(colour = "black")
  )
composition_per_group_C

# plot composition per sample for parkinson disease
composition_per_sample_PD <- plot_composition(ps.top20_PD, 
                                             otu.sort = "abundance",
                                             group_by = "Disease")+
  geom_col(width = 1)+
  scale_fill_viridis_d("Genus")+
  theme_classic()+
  labs(tag = "C")+
  theme(
    axis.text.x = element_text(colour = "black", size = 13, face = "plain", angle = 90, vjust = .4),
    axis.text.y = element_text(colour = "black", size = 13, face = "plain"),
    axis.title.x = element_text(colour = "black", size = 16, face="bold"),
    axis.title.y = element_text(colour = "black", size = 16, face = "bold"),
    strip.text.x = element_text(colour = "black", size = 14, face = "bold"),
    legend.title = element_text(colour = "black", size = 14, face = "plain"),
    legend.text = element_text(colour = "black", size = 10, face = "plain"),
    line = element_line(linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position = "none",
    plot.tag = element_text(hjust=0, size=25, color = "black"),
    axis.ticks = element_line(colour = "black")
  )
composition_per_sample_PD

# plot composition per group for parkinson disease
composition_per_group_PD <- plot_composition(ps.top20_PD, 
                                            otu.sort = "abundance",
                                            average_by = "Disease")+
  geom_col(width = 0.8)+
  scale_fill_viridis_d("Genus")+
  theme_classic()+
  labs(tag = "D")+
  scale_x_discrete(labels= "Average abundance")+
  theme(
    axis.text.x = element_text(colour = "black", size = 13, face = "plain", vjust = -15),
    axis.text.y = element_text(colour = "black", size = 13, face = "plain"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(colour = "black", size = 16, face = "bold"),
    strip.text.x = element_text(colour = "black", size = 14, face = "bold"),
    legend.title = element_text(colour = "black", size = 14, face = "plain"),
    legend.text = element_text(colour = "black", size = 10, face = "plain"),
    line = element_line(linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    plot.tag = element_text(hjust=0, size=25, color = "black"),
    axis.ticks = element_line(colour = "black")
  )
composition_per_group_PD

# 3.4 Visual analysis: Differential abundance -----------------------------

# plot differential abundance to investigate potential shifts driven by host health (i.e., Control or Parkinson disease)
# convert phyloseq object to DESeq2 format
dds = phyloseq_to_deseq2(ps, ~ Disease)
dds = DESeq(dds, test="Wald", fitType="parametric")

# get results 
res = results(dds, cooksCutoff = FALSE)

# set significance threshold and create significance table for each ASV
ALPHA = 0.05

sigtab = res[which(res$padj < ALPHA), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)

# Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)

# create Order column
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)

# create genus column
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
sigtab

# clean up sigtab and add standard error of the mean
sigtab_clean <- sigtab %>%
  group_by(Genus, Order) %>%
  summarise(
    mean_log2FoldChange = mean(log2FoldChange),
    sd_value = sd(log2FoldChange),
    n = n(),
    sem_value = sd_value / sqrt(n)
  )

sigtab_clean

# plot differential abundance 
diff_abundance <- ggplot(sigtab, aes(x = reorder(Genus, log2FoldChange), y = log2FoldChange)) + 
  geom_rect(data=sigtab %>% distinct(Order), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax =0), inherit.aes = F, fill = "cornflowerblue", alpha=0.4)+
  geom_rect(data=sigtab %>% distinct(Order), aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf), inherit.aes = F, fill = "darkorange1", alpha=0.4)+
  geom_point(aes(fill = Order), shape = 21, size = 6, stroke = 1, color = "black", alpha = 0.8) +
  facet_wrap(~ Order, scales = "free_y", ncol = 1) +
  coord_flip()+
  xlab("Genus")+
  #scale_fill_manual(values = c("darkseagreen", "#fde725", "#38598c", "#85d54a", "#2bb07f", "#25858e", "#51c56a", "grey", "#1e9b8a"))+
  scale_fill_manual(values = c("#B39DDB", "#fde725", "#38598c", "#85d54a", "#00BCD4", "brown1", "#51c56a", "grey", "#1e9b8a"))+
  geom_hline(yintercept = 0, color = "black", linewidth=0.75)+
  theme_classic()+
  labs(tag = "A")+
  theme(
    axis.text.x = element_text(colour = "black", size = 13, face = "plain"),
    axis.text.y = element_text(colour = "black", size = 13, face = "plain"),
    axis.title.x = element_text(colour = "black", size = 16, face="bold"),
    axis.title.y = element_text(colour = "black", size = 16, face = "bold"),
    legend.title = element_text(colour = "black", size = 14, face = "plain"),
    legend.text = element_text(colour = "black", size = 10, face = "plain"),
    line = element_line(linewidth = 1),
    panel.spacing = unit(0, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    plot.tag = element_text(hjust=0, size=25, color = "black")
  )


# control bar width across facets
diff_abundance_adj <- diff_abundance + ggforce::facet_col(~ Order, space = "free", scales = "free_y")
diff_abundance_adj

# 3.5 Visual analysis: Beta diversity -------------------------------------

# transform/normalize data and select ordination method and distance measure 
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

# plot beta diversity
beta <- plot_ordination(ps.prop, ord.nmds.bray, color="Disease", title="Bray NMDS")

# get rid of points
beta$layers <- beta$layers[-1]
beta$layers

# plot NMDS to visualize putative differences between Control and Parkinson disease samples
beta_diversity <- beta +
  stat_ellipse(aes(fill = Disease), geom = "polygon", alpha = 0.5, color = "black", linewidth = 0.5) + # add ellipses to better visualize (dis-)similarities
  labs(fill="Group")+ # new legend title
  scale_color_manual(values = c("cornflowerblue", "darkorange1"))+ # add custom colors
  scale_fill_manual(values = c("cornflowerblue", "darkorange1"))+
  geom_point(aes(fill = Disease), shape = 21, color = "black", stroke = 1, size = 4, alpha = 0.6) + # adjust point size, fill and stroke
  guides(color = "none")+ # get rid of one legend
  theme_classic()+
  labs(tag = "B")+
  theme( # adjust looks
    axis.text.x = element_text(colour = "black", size = 13, face = "plain"),
    axis.text.y = element_text(colour = "black", size = 13, face = "plain"),
    axis.title.x = element_text(colour = "black", size = 16, face="bold"),
    axis.title.y = element_text(colour = "black", size = 16, face = "bold"),
    legend.title = element_text(colour = "black", size = 14, face = "plain"),
    legend.text = element_text(colour = "black", size = 10, face = "plain"),
    plot.title = element_blank(),
    line = element_line(linewidth = 1),
    plot.tag = element_text(hjust=0, size=25, color = "black")
  )

beta_diversity

# extract meta data and OTU table, convert to df
ps_meta_df <- data.frame(ps.prop@sam_data)
ps.prop_df <- data.frame(otu_table(ps.prop))

# compute bray-curtis distance matrix
bray_distance <- vegdist(ps.prop_df, method = "bray")

# PERMANOVA
adonis_result <- adonis2(bray_distance ~ Disease, data = ps_meta_df, permutations = 9999)
print(adonis_result)

# check for homogeneity of dispersion between groups
dispersion <- betadisper(bray_distance , ps_meta_df$Disease)
permutest(dispersion)

# 3.6 Visual analysis: Arrange plots --------------------------------------

# Unique genera per sample and group
unique_genera_per_sample_and_group <- (per_sample/per_group)
unique_genera_per_sample_and_group

# save plot 
ggsave("unique_genera.png", units="cm", width=35, height=35, dpi=800)
dev.off()

# Microbiome composition
microbiome_composition <- (composition_per_sample_C+composition_per_group_C+composition_per_sample_PD+composition_per_group_PD)+
  plot_layout(axes = "collect", guides = "keep", widths = c(2, 0.5, 2, 0.5))

microbiome_composition

# save plot 
ggsave("microbiome_composition.png", units="cm", width=40, height=25, dpi=800)
dev.off()

# Differential abundance and beta diversity
differential_abundance_and_beta_diversity <- (diff_abundance_adj+beta_diversity)+plot_layout(guides = "collect", widths = c(1.5, 2))
differential_abundance_and_beta_diversity

# save plot 
ggsave("differential_abundance_and_beta_diversity.png", units="cm", width=40, height=35, dpi=800)
dev.off()


# 4.0 Statistical evaluation ----------------------------------------------

# ALPHA DIVERSITY
# select index and associated values
alpha_data <- alpha$data[6:7] %>%
  group_by(variable)

head(alpha_data)

# loop through data and plot histograms + qq-plots
for (i in unique(alpha_data$variable)){ # get unique index names found in data set
  sub <- subset(alpha_data, variable == i)
  par(mfrow = c(1, 2))
  hist(sub$value, main = i)
  qqnorm(sub$value, main = i)
  qqline(sub$value)
}

# create subsets
chao <- subset(alpha$data, variable == "Chao-1 Index")
simpson <- subset(alpha$data, variable == "Inverse Simpson Index")
shannon <- subset(alpha$data, variable == "Shannon Index")

# students t-test for Chao-1 Index
t_test(value ~ Disease, ref.group = "Control", var.equal = T, paired =F, detailed = T, data = chao, p.adjust.method = "fdr")

# Wilcox-Rank-Sum test for Shannon and inverse Simpson Index
rstatix::wilcox_test(value ~ Disease, ref.group = "Control", exact=T, detailed=T, conf.level = 0.95, data = shannon, p.adjust.method = "fdr")
rstatix::wilcox_test(value ~ Disease, ref.group = "Control", exact=T, detailed=T, conf.level = 0.95, data = simpson, p.adjust.method = "fdr")

# NUMBER OF UNIQUE GENERA PER SAMPLE
# display hist and qq-plot
hist(counts_unique_sample$unique_genera, main = "Unique genera per sample")
qqnorm(counts_unique_sample$unique_genera, main = "Unique genera per sample")
qqline(counts_unique_sample$unique_genera)

# Wilcox-Rank-Sum test 
rstatix::wilcox_test(unique_genera ~ Disease, ref.group = "Control", exact=T, detailed=T, conf.level = 0.95, data = counts_unique_sample, p.adjust.method = "fdr")