# Gut microbiome metagenomics analysis

For this analysis, I have used a subset (n(forward) = 50, n(reverse) = 50) of <a href="https://www.ebi.ac.uk/ena/browser/view/PRJNA494620" >Project: PRJNA494620</a> which provides samples covering a 'Control' group and a group with 'Parkinson disease' with information about the associated gut microbiome. 
I pre-processed the sample primarily using `DADA2` and `cutadapt`. For the subsequent investigation and analysis I used miscellaneous libraries, listed below.

As the folders exceed GitHub's 25MB file size limit, the data can be sent upon request.
## 0.0 Load dependencies
```r
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
```
## 1.0 Set path and import data
```r
# 1.0 Set path and import data -------------------------------------------------------------

# set path 
setwd("insert path")
path <- getwd()
list.files(path)

# use `sort` and `pattern` to ensure that forward and reverse reads are matched correctly
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))  # raw forward reads
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))  # raw reverse reads
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) # define raw sample names and get rid of information about forward or reverse reads
```
### 1.1 Cut primer sequences using `Biostrings`, `DECIPHER` and `cutadapt` 
As primer sequences were not yet removed, the combination of named libraries can help to also remove ambiguous primer sequences, containing e.g., N, V, ..., to improve merging forward and reverse reads later on. 
```r
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
```
**Key output:** 

Investigating the reads.
```r
> head(sread(readFastq(fnFs[1])))
DNAStringSet object of length 6:
    width seq
[1]   250 CCTACGGGGGGCTGCAGTGGGGAATCTTCCGCAATGGACGAAAGTCTGACGGAGCAACGCCGCGTGAGTGATGAAGGA...ACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCGAGCGTTGTCCGGAATTATTGGGCGTAAAGAGCATGTAGGCGGC
[2]   250 CCTACGGGTGGCAGCAGTGGGGAATCTTCCGCAATGGACGAAAGTCTGACGGAGCAACGCCGCGTGAGTGATGAAGGA...ACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCGAGCGTTGTCCGGAATTATTGGGCGTAAAGAGCATGTAGGCGGA
[3]   250 CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGACGCAAGTCTGAACCAGCCATGCCGCGTGCAGGAAGACGGC...CAGCAGCCGCGGTAATACGGAGGATTCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTTTGATA
[4]   250 CCTACGGGTGGCAGCAGTGAGGAATATTGGTCAATGGACGAGAGTCTGAACCAGCCAAGTAGCGTGAAGGATGACTGC...CAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGCGGACGCTTA
[5]   250 CCTACGGGTGGCAGCAGTCGAGAATCATTCACAATGGGGGAAACCCTGATGGTGCGACGCCGCGTGGGGGAATGAAGG...ATACAGAGGTCTCAAGCGTTGTTCGGAATCACTGGGCGTAAAGCGTGCGTAGGCTGTTTCGTAAGTCGTGTGTGAAA
[6]   250 CCTACGGGAGGCAGCAGTGAGGAATATTGGTCAATGGACGCAAGTCTGAACCAGCCATGCCGCGTGCAGGATGACGGC...CAGCAGCCGCGGTAATACGGAGGATTCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTTAGATA
```
Find the number of primer fitting the given sequence.
```r
> cat("Found", length(Fs_chr), "forward primer sequences:", Fs_chr)
Found 8 forward primer sequences: CCTACGGGAGGCAGCAG CCTACGGGCGGCTGCAG CCTACGGGGGGCAGCAG CCTACGGGTGGCTGCAG CCTACGGGAGGCTGCAG CCTACGGGCGGCAGCAG CCTACGGGGGGCTGCAG CCTACGGGTGGCAGCAG
> cat("Found", length(Rs_chr), "reverse primer sequences:", Rs_chr)
Found 9 reverse primer sequences: GACTACAAGGGTATCTAATCC GACTACTCGGGTATCTAATCC GACTACCGGGGTATCTAATCC GACTACACGGGTATCTAATCC GACTACTGGGGTATCTAATCC GACTACCAGGGTATCTAATCC GACTACAGGGGTATCTAATCC GACTACTAGGGTATCTAATCC GACTACCCGGGTATCTAATCC
```
Some output from `cutadapt`
```r
=== Summary ===

Total read pairs processed:            424,236
  Read 1 with adapter:                 408,271 (96.2%)
  Read 2 with adapter:                 415,472 (97.9%)

== Read fate breakdown ==
Pairs discarded as untrimmed:           23,893 (5.6%)
Pairs written (passing filters):       400,343 (94.4%)

Total basepairs processed:   212,118,000 bp
  Read 1:   106,059,000 bp
  Read 2:   106,059,000 bp
Total written (filtered):    184,702,568 bp (87.1%)
  Read 1:    93,279,919 bp
  Read 2:    91,422,649 bp
```
### 1.2 Visualize quality profiles, `filterAndTrim()` of cut files
After cutting out the primer sequences from each sample, the files read quality was assessed to determine if or how much of the reads must be trimmed of. Here, I used a custom function to loop through each provided sample to
- visualize each sample and annotate the plot,
- cut the files based on a provided threshold of mean read quality that, if undershot, saves the x-axis information (i.e., Cycle) in a separat data frame to be fed to the `DADA2` function `filterAndTrim()`.

The files were filtered and trimmed, the amplicon length was calculated to report the overlap between forward and reverse reads. Moreover, the percentage of filtered reads was also printed out to get a better idea of the filtering process and how much data was lost.
```r
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
```
**Key output:**  
`where_to_trim_Fs()`

No trimming needed for SRR16873586 based on set threshold of 27 .  

<img width="819" height="588" alt="image" src="https://github.com/user-attachments/assets/e8e15298-0292-4ab8-b814-27a3503b918d" />  

`where_to_trim_Rs()`  

Suggested trimming of SRR16873586 based on set threshold of 26 at 192 cycles with a value of 25.64 .

<img width="819" height="588" alt="image" src="https://github.com/user-attachments/assets/9404c38e-1663-40cf-a0b9-1c5875890f89" />


Amplicon length to calculate the overlap between forward and reverse reads and resulting data frame `trim_at_df_flex` based on the output of the function mentioned above.

```r
> amplicon_length <- calc_amplicon_length(forward_primer_binding_site="341", reverse_primer_binding_site = "805", Fs_primer="CCTACGGGNGGCWGCAG", Rs_primer="GACTACHVGGGTATCTAATCC")
The amplicon length is 427 bp.
```
<img width="817" height="586" alt="image" src="https://github.com/user-attachments/assets/8e85e380-3179-4e31-bd81-1e6cebc1917d" />

Resulting data frame with reads in and reads out as well as the percentage of reads after filtering and trimming.
```r
> out_df
            reads.in reads.out reads.out (%)
SRR16873583   946589    688667      72.75248
SRR16873586   359244    257443      71.66244
SRR16873589   451001    317828      70.47168
SRR16873591  1119159    843238      75.34568
SRR16873592   369252    266119      72.06975
SRR16873595   313460    224377      71.58074
SRR16873597   370881    261357      70.46923
SRR16873608   306961    210541      68.58884
SRR16873609   334484    235620      70.44283
SRR16873620   359588    229279      63.76158
SRR16873624   441271    312894      70.90745
SRR16873627   325551    201633      61.93592
SRR16873630   319962    218924      68.42188
SRR16873632   351056    228938      65.21410
SRR16873633   390921    280972      71.87437
SRR16873635  1003739    647530      64.51179
SRR16873649   278707    195624      70.18984
SRR16873650   309042    218516      70.70754
SRR16873653   733604    466422      63.57953
SRR16873659   331999    232147      69.92401
SRR16873660   640021    418145      65.33301
SRR16873662   428937    285483      66.55593
SRR16873663   349058    252975      72.47363
SRR16873668  1055396    703015      66.61149
SRR16873675   298656    211969      70.97430
SRR16873679   292117    204839      70.12225
SRR16873680   790492    560883      70.95366
SRR16873682   558432    387950      69.47131
SRR16873683   702476    460151      65.50416
SRR16873687   472139    343330      72.71799
SRR16873692   419216    295685      70.53285
SRR16873695   479462    346939      72.36006
SRR16873696   528144    376921      71.36709
SRR16873697   355750    235015      66.06184
SRR16873700   692737    449048      64.82229
SRR16873706   360815    261852      72.57237
SRR16873708   718338    463687      64.54998
SRR16873709   371608    269313      72.47234
SRR16873710   303692    214667      70.68576
SRR16873713  1221601    813109      66.56093
SRR16873746   700259    468638      66.92352
SRR16873747   587470    373395      63.55984
SRR16873751   529416    384808      72.68537
SRR16873756   305372    214231      70.15411
SRR16873761   502831    336931      67.00681
SRR16873765   310231    215519      69.47049
SRR16873767   674841    450723      66.78951
SRR16873768   875231    574882      65.68346
SRR16873773   430300    285038      66.24169
SRR16873779   400343    261007      65.19584
```

### 1.3 Learn errors, merge forward and reverse reads and remove chimeras

```r
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
```
**Key output:**
```r
> # remove chimeras
> seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE) 
Identified 34198 bimeras out of 39980 input sequences.
> dim(seqtab.nochim) # number of ASVs left after getting rid of chimeras
[1]   50 5782
> 
> sum(seqtab.nochim)/sum(seqtab) # 0-1; proportion of reads after chimera removal.
[1] 0.9504996
```
### 1.4 Pipeline summary
This section serves as a summary of the efforts to pre-process, trim, filter, and clean up the loaded data. The resulting data frame represents
- absolute input reads,
- absolute filtered reads,
- filtered reads (%),
- denoised reads,
- merged reads,
- non-chimera reads,
- merge rate (%).

A merge rate of 75% or higher is considered sufficient and can be used with confidence. 
```r
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
```
<img width="817" height="586" alt="image" src="https://github.com/user-attachments/assets/749c3bbe-bc0a-47e9-afc3-c0b7a93564f3" />

### 1.5 Handle samples with low merge rate
To exclude samples with a low merge rate and thus a sample that would yield results with a low confidence, the data needs to be filtered. 
```r
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
```
**Key output:**  
Seven samples were excluded.
```r
+   print(nonchim)
+   length(nonchim)
+ }
SRR16873583 SRR16873589 SRR16873595 SRR16873608 SRR16873609 SRR16873620 SRR16873624 SRR16873627 SRR16873630 SRR16873632 
     624674      292752      196093      201728      221573      197069      269804      187935      180084      199945 
SRR16873633 SRR16873635 SRR16873650 SRR16873653 SRR16873659 SRR16873660 SRR16873662 SRR16873663 SRR16873668 SRR16873675 
     238638      574396      157264      419902      193017      379888      275267      239491      593216      190621 
SRR16873679 SRR16873680 SRR16873682 SRR16873683 SRR16873695 SRR16873696 SRR16873697 SRR16873700 SRR16873706 SRR16873708 
     189346      425137      328316      428987      328417      350109      226293      423864      224670      400211 
SRR16873709 SRR16873710 SRR16873713 SRR16873746 SRR16873747 SRR16873751 SRR16873756 SRR16873761 SRR16873765 SRR16873767 
     244549      203552      767593      412110      345151      323414      187240      308561      190704      393844 
SRR16873768 SRR16873773 SRR16873779 
     505401      266326      240894 
[1] 43
```
### 1.6 Assign taxonomy
```r
# 1.6 Assign taxonomy -----------------------------------------------------

# assign taxonomy(Genus and/or Species)
taxa <- assignTaxonomy(seqtab.nochim_cleaned, "./silva_nr99_v138.2_toSpecies_trainset.fa.gz", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```
**Key output:**
```r
> head(taxa.print)
     Kingdom    Phylum              Class              Order                Family               Genus                  
[1,] "Bacteria" "Bacteroidota"      "Bacteroidia"      "Bacteroidales"      "Bacteroidaceae"     "Bacteroides"          
[2,] "Bacteria" "Bacteroidota"      "Bacteroidia"      "Bacteroidales"      "Bacteroidaceae"     "Bacteroides"          
[3,] "Bacteria" "Bacteroidota"      "Bacteroidia"      "Bacteroidales"      "Bacteroidaceae"     "Bacteroides"          
[4,] "Bacteria" "Verrucomicrobiota" "Verrucomicrobiia" "Verrucomicrobiales" "Akkermansiaceae"    "Akkermansia"          
[5,] "Bacteria" "Bacteroidota"      "Bacteroidia"      "Bacteroidales"      "Bacteroidaceae"     "Bacteroides"          
[6,] "Bacteria" "Bacillota"         "Negativicutes"    "Acidaminococcales"  "Acidaminococcaceae" "Phascolarctobacterium"
     Species      
[1,] "vulgatus"   
[2,] "dorei"      
[3,] "uniformis"  
[4,] "muciniphila"
[5,] "stercoris"  
[6,] "faecium"
```
## 2.0 Import meta data 
To get insights, if or how the microbiome changed depending on being diagnosted with Parkinson disease, the meta data needs to be imported. A for loop can loop through provided .xml files and append the raw values to an empty data frame to collect all the important information needed for the subsequent analysis.
```r
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
```
**Key output:**
```r
                      Disease    Sex Tissue Age
SRR16873583 Parkinson disease   male  stool  59
SRR16873586 Parkinson disease female  stool  66
SRR16873589 Parkinson disease   male  stool  66
SRR16873591 Parkinson disease   male  stool  70
SRR16873592 Parkinson disease female  stool  74
SRR16873595 Parkinson disease   male  stool  67
SRR16873597 Parkinson disease   male  stool  69
SRR16873608 Parkinson disease female  stool  72
SRR16873609 Parkinson disease   male  stool  75
SRR16873620 Parkinson disease female  stool  77
SRR16873624 Parkinson disease female  stool  59
SRR16873627 Parkinson disease female  stool  62
SRR16873630 Parkinson disease   male  stool  64
SRR16873632 Parkinson disease female  stool  62
SRR16873633 Parkinson disease   male  stool  48
SRR16873635 Parkinson disease   male  stool  69
SRR16873649 Parkinson disease   male  stool  59
SRR16873650 Parkinson disease female  stool  55
SRR16873653 Parkinson disease female  stool  58
SRR16873659 Parkinson disease female  stool  53
SRR16873660 Parkinson disease female  stool  77
SRR16873662 Parkinson disease   male  stool  60
SRR16873663 Parkinson disease   male  stool  66
SRR16873668 Parkinson disease   male  stool  64
SRR16873675 Parkinson disease   male  stool  62
SRR16873679 Parkinson disease female  stool  46
SRR16873680 Parkinson disease   male  stool  68
SRR16873682           Control   male  stool  71
SRR16873683           Control   male  stool  70
SRR16873687           Control female  stool  58
SRR16873692           Control female  stool  64
SRR16873695           Control female  stool  63
SRR16873696           Control female  stool  52
SRR16873697           Control female  stool  74
SRR16873700           Control female  stool  66
SRR16873706           Control female  stool  67
SRR16873708           Control   male  stool  59
SRR16873709           Control   male  stool  61
SRR16873710           Control   male  stool  64
SRR16873713           Control   male  stool  72
SRR16873746           Control female  stool  67
SRR16873747           Control female  stool  63
SRR16873751           Control   male  stool  63
SRR16873756           Control   male  stool  67
SRR16873761           Control female  stool  56
SRR16873765           Control   male  stool  50
SRR16873767           Control female  stool  59
SRR16873768           Control female  stool  62
SRR16873773           Control   male  stool  48
SRR16873779           Control female  stool  65
```
## 3.0 Prepare data and data wrangling
The created files holding information about abundances, sequences, ... get combined in a single `phyloseq` object named `ps`. Moreover, the file can be saved and loaded if desired. The data can then be sorted and selected for the 20 most abundant taxa (e.g., genus) for relative abundance plots.
```r
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
```
**Key output:**  

Top 20 most abundant genera in Control samples
```r
> ps.top20_C
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 12 taxa and 21 samples ]
sample_data() Sample Data:       [ 21 samples by 4 sample variables ]
tax_table()   Taxonomy Table:    [ 12 taxa by 7 taxonomic ranks ]
```
Top 20 most abundant genera in Parkinson disease samples
```r
> ps.top20_PD
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 12 taxa and 22 samples ]
sample_data() Sample Data:       [ 22 samples by 4 sample variables ]
tax_table()   Taxonomy Table:    [ 12 taxa by 7 taxonomic ranks ]
```
### 3.1 Visual analysis: Alpha diversity 
```r
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
```
The Chao-1 displayed that both groups exhibited a similar alpha diversity. However, the control presented a higher variance and median compared to Parkinson disease samples. There was no statistically significant difference found between tested groups (students t-test, p > 0.05). Therefore, all samples showed a similar species richness, with no difference between patients.

The Inverse Simpson index viewed that control and Parkinson disease samples showed a similar median value, but control samples displayed a higher variance and thus a greater diversity. The Wilcoxon rank sum test revealed no statistically significant difference between the groups (p > 0.05).  

The Shannon index measures illustrated that control samples returned a greater variance with a lower alpha diversity in some samples compared to Parkinson disease patients. Diseased patients displayed a lower median value but also a lower variance, suggesting a more even species distribution (i.e., richness and eveness) across samples. No statistically significant difference was found between sampled groups (Wilcoxon rank sum test, p > 0.05).
<img width="11023" height="7874" alt="alpha_diversity" src="https://github.com/user-attachments/assets/bf8161bf-74d9-4873-8497-a3a1c7b51eb3" />
<br />
<br />
### 3.2 Visual analysis: Unique genera per sample and group
To further investigate the taxon diversity, the number of unique genera per sample and group was visually analyzed. However, the number of unique genera per sample does not necessarily reflect the total number. Thus, the total number of unique genera was visually analyzed as well with an additional venn diagram demonstrating the share of identical and unique genera with an abundance > 0.
```r
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
Parkinson_genera <- ps_df_selected$`Parkinson disease`

# unique genera found in control
unique_genera_control <- setdiff(control_genera, Parkinson_genera)
unique_genera_control

# unique genera found in Parkinson disease samples
unique_genera_pd <- setdiff(Parkinson_genera, control_genera)
unique_genera_pd

# plot venn diagram
venn_plot <- ggvenn::ggvenn(list(
  Control = control_genera,
  Parkinson = Parkinson_genera),
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
```
As displayed in the alpha diversity measures, the control exhibited more samples with a lower diversity than Parkinson disease patients (A). However, the mean and median value between groups is very similar and thus, no statistically significant difference could be found (Wilcoxon rank sum test, p > 0.05)(B). 

The number of unique genera per group (C) showed that the control group contains 262 genera while the Parkinson disease group contains only 241 unique members. However, as these single values do not allow for a statistical evaluation, this insight may only suggest, that there is a greater variance in the control group. The venn diagram (D) supports this, as both groups share 216 genera (75%) and differ moderately from one another. 
<img width="11023" height="11023" alt="unique_genera" src="https://github.com/user-attachments/assets/2ebaea48-d360-4557-925e-9a0585124480" />
<br />
<br />
### 3.3 Visual analysis: Microbiome composition
To investigate putative differences in abundance between samples and groups, visualizing the microbiome composition can be of help. Here, the top 20 most abundant genera were used to visually display the fraction of each genus per sample and the average abundance within each group. 
```r
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

# plot composition per sample for Parkinson disease
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

# plot composition per group for Parkinson disease
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
```
(A) The control group displayed differences between samples, especially in genera with a lower overall abundance (e.g., Acidaminococcus, Parabacteroides). However, many samples showed a high fraction of Bacteroides, Megasphaera and Alistipes, making them the most dominant genus. 
(B) The average abundance across samples illustrated that Bacteroides make up almost 50% of the total abundance, followed by Megasphaera, Alistipes and Phascolarctobacterium. The less abundant genera displayed a similar average fraction across samples.  

(C) Patients with Parkinson disease displayed a high similarity to the control group but showed a higher fraction of moderately abundant genera (e.g., Megasphaera, Akkermansia). 
(D) As presented in control samples, Bacteroides is the most abundant genus, followed by Megasphaera. However, other genera occur at different abundance measures with higher fractions in Megasphaera and Akkermansia, thus ordered differently in the legend. 
Overall, Parkinson patients represented a similar fraction of Bacteroides and Megasphaera but seem to yield different fractions and abundances of the other genera, suggesting a shift in the microbiome composition. 
<img width="12598" height="7874" alt="microbiome_composition" src="https://github.com/user-attachments/assets/429539fa-2d90-4fa8-9cfc-f253bf63021a" />
<br />
<br />
### 3.4 Visual analysis: Differential abundance and beta diversity

```r
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
```
**Key output:**
```r
> # PERMANOVA
> adonis_result <- adonis2(bray_distance ~ Disease, data = ps_meta_df, permutations = 9999)
> print(adonis_result)
Permutation test for adonis under reduced model
Permutation: free
Number of permutations: 9999

adonis2(formula = bray_distance ~ Disease, data = ps_meta_df, permutations = 9999)
         Df SumOfSqs      R2      F Pr(>F)
Model     1   0.3512 0.02395 1.0062 0.4443
Residual 41  14.3120 0.97605              
Total    42  14.6632 1.00000
```
```r
> # check for homogeneity of dispersion between groups
> dispersion <- betadisper(bray_distance , ps_meta_df$Disease)
> permutest(dispersion)

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df   Sum Sq    Mean Sq      F N.Perm Pr(>F)
Groups     1 0.000181 0.00018122 0.0691    999  0.805
Residuals 41 0.107541 0.00262294
```
(A) Data points with positive log2FoldChange values represent an enrichment of the associated genus in Parkinson disease patients, represented by an orange background. Data points with a negative log2FoldChange value represent an enrichment of the associated genus in the control group, represented by a blue background.  

Parkinson disease patients displayed a statistically significant enrichment in Holdemanella of the order Erysipelotrichales, unknown genera as well as most of the other genera of the order Bacteroidales. However, there was also a statistically signigicant enrichment found for Bacteroides, Parabacteroides and Odoribacter in the control group. The genus of Phascolarctobacterium and many genera of the order Oscillospirales depicted a vast increase in abundance in the Parkinson disease group. On the other hand, unknown genera and Intestinomonas and Gemmiger were enriched in the control group. In addition, unknown genera and Lachnospiraceae UCG-010 rose in abundance compared to control. In control, Defluviitaleaceae UCG-011 and Agathobacter gained in abundance. Megasphaera and Akkermansia displayed an enrichment in Parkinson disease patients while Christensenellaceae R-7 group and Bifidobacterium exhibited increments in abundance within the control group. 

Overall, for both groups statistically significant gains in abundance were found, suggesting a shift in microbiome composition when diagnosed with Parkinson disease. 

(B) The beta-diversity visualization, using the bray-curtis distance matrix, illustrating the difference between control and Parkinson disease samples. However, no statistically significant differences were found as there is a large overlap between both groups, represented by ellipses. Thus, suggesting no difference in microbial composition between mentioned groups. 
<img width="12598" height="11023" alt="differential_abundance_and_beta_diversity" src="https://github.com/user-attachments/assets/d25aeb57-3bd9-4d24-ae4b-6e35ace88ade" />
<br />
<br />

## 4.0 Statistical analysis

```r
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
hist(phylum_counts_unique_sample$unique_genera)
qqnorm(phylum_counts_unique_sample$unique_genera)
qqline(phylum_counts_unique_sample$unique_genera)

# Wilcox-Rank-Sum test 
rstatix::wilcox_test(unique_genera ~ Disease, ref.group = "Control", exact=T, detailed=T, conf.level = 0.95, data = complete, p.adjust.method = "fdr")
```
**Key output:**  

Alpha-diversity 
<img width="984" height="715" alt="image" src="https://github.com/user-attachments/assets/e9ab4c7b-1708-4c6c-be42-8a8600cb7377" />
<img width="984" height="715" alt="image" src="https://github.com/user-attachments/assets/a70967fe-9cf5-4355-a786-845788c38ee5" />
<img width="984" height="715" alt="image" src="https://github.com/user-attachments/assets/861ca68c-b9ac-4f1a-8c65-18cbccaf53e1" />

Unique genera per sample
<img width="984" height="715" alt="image" src="https://github.com/user-attachments/assets/42b67174-611c-474d-9a0e-10645e1b2998" />

## 5.0 Conclusion
There was a similar alpha-diversity found for both tested groups with a larger variance for control samples. Furthermore, a similar number of unique genera was displayed with a high overlap between control and Parkinson disease with only a moderate difference between groups, although based on no statistical evaluation. 

However, the investigation of the differential abundance showed, that with the diagnosis of Parkinson disease, a shift in the enrichment of various genera could be observed. Based on beta-diversity measures, no difference in the microbial compositions could be discovered.  


