# Amplicon metagenomics

## Process raw data using DADA2

The following code was used to process raw sequencing data (`*.fastq` files) to generate ASV tables. Metadata and the output `ASV_table.rds` are included in the `data` folder, and can be used to replicate the analyses below.

```R
library("tidyverse")
library("gridExtra")
library("dada2")
library("phyloseq")
library("DECIPHER")
library("phangorn")
library("ape")
library("seqinr")

n_cores <- 16

dir.create('SILVA', showWarnings = FALSE, recursive = TRUE)
download.file("https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz", file.path('SILVA', 'silva.fa.gz'))

dir.create('PR2', showWarnings = FALSE, recursive = TRUE)
download.file("https://github.com/pr2database/pr2database/releases/download/v4.14.0/pr2_version_4.14.0_SSU_dada2.fasta.gz", file.path('PR2', 'pr2.fa.gz'))

process.files <- function(indir, filter_dir, outdir, tax, tax_key, metadata_file, tree, asv.final.table){
 
  dir.create(filter_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  fnFs <- sort(list.files(indir, pattern="_1.fq.gz", full.names = TRUE))
  fnRs <- sort(list.files(indir, pattern="_2.fq.gz", full.names = TRUE))
  sampleIDs <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  
  filtFs <- file.path(filter_dir, paste0(sampleIDs, "_F_filt.fastq.gz"))
  filtRs <- file.path(filter_dir, paste0(sampleIDs, "_R_filt.fastq.gz"))
  names(filtFs) <- sampleIDs
  names(filtRs) <- sampleIDs
  
  cat("Filterning and trimming ...", format(Sys.time(), "%c"), "\n")
  
  filter_results <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                                  rm.phix = FALSE,
                                  multithread = n_cores, 
                                  compress = FALSE, verbose = TRUE)
  
  
  exists <- file.exists(fnFs) & file.exists(filtFs)
  filtFs <- filtFs[exists]
  exists <- file.exists(fnRs) & file.exists(filtRs)
  filtRs <- filtRs[exists]
  
  cat("Learning error rates ...", format(Sys.time(), "%c"), "\n")
  
  errF <- learnErrors(filtFs, multithread = n_cores, verbose = TRUE)
  errR <- learnErrors(filtRs, multithread = n_cores, verbose = TRUE)
  
  cat("Dereplicating sequences ...", format(Sys.time(), "%c"), "\n")
  
  dadaFs <- dada(filtFs, err=errF, pool = FALSE, multithread = n_cores)
  dadaRs <- dada(filtRs, err=errR, pool = FALSE, multithread = n_cores)
  
  cat("Merging sequences ...", format(Sys.time(), "%c"), "\n")
  
  mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
  
  cat("Creating ASV table ...", format(Sys.time(), "%c"), "\n")
  
  seqtab_all <- makeSequenceTable(mergers)
  
  cat("Looking for chimeras ...", format(Sys.time(), "%c"), "\n")
  
  seqtab <- removeBimeraDenovo(seqtab_all,
                               method = 'consensus',
                               multithread = n_cores,
                               verbose = TRUE)
  
  cat("Assigning taxonomy ...", format(Sys.time(), "%c"), "\n")
  
  taxa <- assignTaxonomy(seqtab, tax_key, multithread = n_cores)
  
  otutable <- seqtab
  colnames(otutable) <- paste('ASV', 1:ncol(seqtab), sep = '_')
  row.names(taxa) <- paste('ASV', 1:ncol(seqtab), sep = '_')
  
  asv_seqs <- colnames(seqtab)
  asv_headers <- paste('ASV', 1:ncol(seqtab), sep = '_')
  asv_fasta <- c(rbind(asv_headers, asv_seqs))
  
  cat("Aligning sequences ...", format(Sys.time(), "%c"), "\n")
  
  seqs <- getSequences(seqtab)
  names(seqs) <- asv_headers 
  alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA, processors = n_cores)
  phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
  
  cat("Calculating tree ...", format(Sys.time(), "%c"), "\n")
  
  dna_dist <- dist.ml(phang.align, model="JC69")
  asv_UPGMA <- upgma(dna_dist)
  fit <- pml(asv_UPGMA, phang.align)
  fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")
  tree <- bootstrap.pml(fitJC, bs=1, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))
  
  cat("Creating phyloseq object ...", format(Sys.time(), "%c"), "\n")
  
  ps <- phyloseq(otu_table(otutable, taxa_are_rows = FALSE),
                 tax_table(taxa),
                 tree[[1]])
  
  return(ps)
}


ps.16s <- process.files(
   indir <- '16s/data',
   filter_dir <- '16s/fastq_filtered',
   outdir <- '16s/ASVs',
   taxdir <- 'SILVA',
   tax_key <- 'SILVA/silva.fa.gz',
   metadata_file <- 'metadata_16s.txt',
   asv.final.table <- 'ASV_table_16s.rds',
   tree <- 'tree.16s.tre'
)

ps.18s <- process.files(
  indir <- '18s/data',
  filter_dir <- '18s/fastq_filtered',
  outdir <- '18s/ASVs',
  taxdir <- 'PR2',
  tax_key <- 'PR2/pr2.fa.gz',
  metadata_file <- 'metadata_18s.txt',
  asv.final.table <- 'ASV_table_18s.rds',
  tree <- 'tree.18s.tre'
)

rm(list=setdiff(ls(), c("ps.16s", "ps.18s")))
save.image(file = 'ASV_table.rds')
```

# Data analysis

## Figure 5b

### Load libraries

```R
library("phyloseq")
library("ggplot2")
library("DESeq2")
library("vegan")
library("ape")
library("plyr") 
library("ggpmisc")
library("dplyr")
library("broom")
library("picante")
library("Rmisc")
library("emmeans")
library("car")
library("lme4")
library("tibble")
library("data.table")
library("limma")
library("microbiome")
library("RColorBrewer")
library("RVAideMemoire")
library("data.table")
library("decontam")
library("ggpubr")
library("cowplot")
```

### Load data

```R
load(file = 'data/ASV_table.rds')
metadata <- read.table("data/metadata_16s.txt", sep = "\t", header = T, row.names = 1)
sample_data(ps.16s) <- metadata
metadata <- read.table("data/metadata_18s.txt", sep = "\t", header = T, row.names = 1)
sample_data(ps.18s) <- metadata
```

### Remove contaminants and cleanup the dataset

```R
remove.cont <- function(GM){
  sample_data(GM)$is.neg <- sample_data(GM)$genotype == "nc"
  contamdf.prev <- isContaminant(GM, method="prevalence", neg="is.neg", threshold = 0.05)
  cont.remove <- subset(contamdf.prev, contaminant == "TRUE")
  cont.remove <- row.names(cont.remove)
  allTaxa = taxa_names(GM)
  allTaxa <- allTaxa[!(allTaxa %in% cont.remove)]
  GM <-  prune_taxa(allTaxa, GM)
  GM <- subset_samples(GM, genotype != "nc")
  return(GM)
}

ps.16s <- remove.cont(ps.16s)
ps.16s <- subset_taxa(ps.16s, Order !="Chloroplast")
ps.16s <- subset_taxa(ps.16s, Family !="Mitochondria")
ps.16s <- filter_taxa(ps.16s, function (x) {sum(x > 0) > 1}, prune=TRUE)

ps.18s <- remove.cont(ps.18s)
ps.18s <- subset_taxa(ps.18s, Class !="Streptophyta")
ps.18s <- filter_taxa(ps.18s, function (x) {sum(x > 0) > 1}, prune=TRUE)
```

### Figure 5b (Sp21) - plot

```R
GM <- subset_samples(ps.16s, genotype == "Sp21")
dist.mat <- phyloseq::distance(GM, method = "bray")
cap_ord <- ordinate(physeq = GM, method = "CAP", distance = dist.mat, formula = ~ cage)
cap_plot1 <- plot_ordination(physeq = GM, ordination = cap_ord, axes = c(1,2)) +
  theme_bw(base_size = 15) +
  stat_ellipse(mapping = aes(color = cage),
               alpha = 0.4,
               type = "norm",
               show.legend=F) +
  geom_point(mapping = aes(color = cage), size = 5) +
  theme(legend.title= element_blank(), 
        legend.background = element_rect(color = NA),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        panel.grid = element_blank()) +
  scale_color_manual(name = "Legend", values=c("#4daf4a", "#e41a1c"), labels = c("Control", "Herbivory"), breaks = c("ctrl", "herb"))
cap_plot1
```

### Figure 5b (Sp21) - PERMANOVA

```R
GM <- subset_samples(ps.16s, genotype == "Sp21")
sampledf <- data.frame(sample_data(GM))
dist.mat <- phyloseq::distance(GM, method = "bray")
perm <- how(nperm = 999)
setBlocks(perm) <- with(sampledf, pond)
set.seed(100)
pmv <- adonis2(dist.mat ~ cage, data = sampledf, permutations = perm)
pmv
```

### Figure 5c (Sp65) - plot

```R
GM <- subset_samples(ps.16s, genotype == "Sp65")
dist.mat <- phyloseq::distance(GM, method = "bray")
cap_ord <- ordinate(physeq = GM, method = "CAP", distance = dist.mat, formula = ~ cage)
cap_plot1 <- plot_ordination(physeq = GM, ordination = cap_ord, axes = c(1,2)) +
  theme_bw(base_size = 15) +
  stat_ellipse(mapping = aes(color = cage),
               alpha = 0.4,
               type = "norm",
               show.legend=F) +
  geom_point(mapping = aes(color = cage), size = 5) +
  theme(legend.title= element_blank(), 
        legend.background = element_rect(color = NA),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        panel.grid = element_blank()) +
  scale_color_manual(name = "Legend", values=c("#4daf4a", "#e41a1c"), labels = c("Control", "Herbivory"), breaks = c("ctrl", "herb"))
cap_plot1
```

### Figure 5c (Sp65) - PERMANOVA

```R
GM <- subset_samples(ps.16s, genotype == "Sp65")
sampledf <- data.frame(sample_data(GM))
dist.mat <- phyloseq::distance(GM, method = "bray")
perm <- how(nperm = 999)
setBlocks(perm) <- with(sampledf, pond)
set.seed(100)
pmv <- adonis2(dist.mat ~ cage, data = sampledf, permutations = perm)
pmv
```

## Table S8

First we define the two function for creating the table. `up` is used to identify ASVs more abundant under herbivory, `down` is used to identify ASVs more abundant under control conditions.

```R
cal.diff.taxa <- function(object, direction){
  diagdds <- phyloseq_to_deseq2(object, ~ 1)
  ts <- counts(diagdds)
  geoMeans = apply(ts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  diagdds = estimateSizeFactors(diagdds, geoMeans=geoMeans)
  diagdds = estimateDispersions(diagdds)
  diagdds$group <- factor(paste0(diagdds$cage))
  design(diagdds) <- ~ group
  dds <-DESeq(diagdds, betaPrior=FALSE, parallel = T)
  c1 <- results(dds, contrast=c("group", "herb", "ctrl"), parallel = T)
  c1 <- as.data.frame(c1)
  c1 <- setDT(c1, keep.rownames = TRUE)[]
  c1 <- c1[,c("rn", "log2FoldChange", "padj")]
  x <- if(direction == "up"){c1 %>% dplyr::filter(log2FoldChange > 0 & padj < 0.05)} else if(direction == "down"){c1 %>% dplyr::filter(log2FoldChange < 0 & padj < 0.05)}
  tax.table <- as.data.frame(tax_table(object))
  tx <- tax.table[which(row.names(tax.table) %in% x$rn),]
  return(tx)
}

create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel'),
                               lengthMenu = list(c(10,25,50,-1),
                                                 c(10,25,50,"All"))))
}
```

ASVs more abundant under herbivory in genotype **Sp21**

```R
GM <- subset_samples(ps.16s, genotype == "Sp21")
cal.diff.taxa(GM, "up") %>% create_dt
```

ASVs more abundant under control conditions in genotype **Sp21**

```R
GM <- subset_samples(ps.16s, genotype == "Sp21")
cal.diff.taxa(GM, "down") %>% create_dt
```

ASVs more abundant under herbivory in genotype **Sp65**

```R
GM <- subset_samples(ps.16s, genotype == "Sp65")
cal.diff.taxa(GM, "up") %>% create_dt
```

ASVs more abundant under control conditions in genotype **Sp65**

```R
GM <- subset_samples(ps.16s, genotype == "Sp65")
cal.diff.taxa(GM, "down") %>% create_dt
```

# PERMANOVA and differential ASVs for 18S

## Sp21
```R
GM <- subset_samples(ps.18s, genotype == "Sp21")
sampledf <- data.frame(sample_data(GM))
dist.mat <- phyloseq::distance(GM, method = "bray")
perm <- how(nperm = 999)
setBlocks(perm) <- with(sampledf, pond)
set.seed(100)
pmv <- adonis2(dist.mat ~ cage, data = sampledf, permutations = perm)
pmv
```

## Sp65
```R
GM <- subset_samples(ps.18s, genotype == "Sp65")
sampledf <- data.frame(sample_data(GM))
dist.mat <- phyloseq::distance(GM, method = "bray")
perm <- how(nperm = 999)
setBlocks(perm) <- with(sampledf, pond)
set.seed(100)
pmv <- adonis2(dist.mat ~ cage, data = sampledf, permutations = perm)
pmv
```
