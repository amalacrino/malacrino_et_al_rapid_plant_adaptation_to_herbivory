# Pool sequencing

## Trim Galore

```Bash
DATADIR=path_to_data
OUTDIR=path_to_output

cd $DATADIR

find -name "*_R1.fq.gz" | rev | cut -c 10- | rev | parallel -j 32 trim_galore --illumina --paired --fastqc --gzip -o $OUTDIR/ {}\_R1.fq.gz {}\_R2.fq.gz
```

## Bowtie

We use Bowtie to map reads to the *S. polyrhiza* reference genome. `BAM` files with reads mapping the genome will be used to estimate genotype frequency (see HAFpipe below), while reads not mapping the references are converted back to `fastq` files and used for inferring the plant-associated microbial community (see Kraken2 below).

```Bash
DATADIR=path_to_TrimGalore_output
REFDIR=path_to_reference_genome
OUTDIR=path_to_output_directory
OUTPUT_mapped=$OUTDIR/mapped
OUTPUT_unmapped=$OUTDIR/host_removed

NTHREADS=32

cd $REFDIR

bowtie2-build --threads $NTHREADS SP_combined.fasta duckweed_reference

cd $DATADIR

for file in *_R1_val_1.fq.gz
do
SAMPLE=$(basename $file "_R1_val_1.fq.gz")
echo "$SAMPLE Bowtie2"
bowtie2 -p $NTHREADS -x $REFDIR/duckweed_reference -1 ${SAMPLE}_R1_val_1.fq.gz -2 ${SAMPLE}_R2_val_2.fq.gz -S $OUTDIR/${SAMPLE}.sam

echo "$SAMPLE conversion from sam to bam"
samtools view -@ $NTHREADS -bS $OUTDIR/${SAMPLE}.sam > $OUTDIR/${SAMPLE}.bam

echo "$SAMPLE select unmapped"
samtools view -@ $NTHREADS -b -f 12 -F 256 $OUTDIR/${SAMPLE}.bam > $OUTDIR/${SAMPLE}_unmapped.bam

echo "$SAMPLE select mapped"
samtools view -@ $NTHREADS -b -f 3 $OUTDIR/${SAMPLE}.bam > $OUTDIR/${SAMPLE}_mapped.bam

echo "$SAMPLE sorting unmapped"
samtools sort -@ $NTHREADS -n $OUTDIR/${SAMPLE}_unmapped.bam -o $OUTDIR/${SAMPLE}_unmapped_sorted.bam

echo "$SAMPLE sorting mapped"
samtools sort -@ $NTHREADS $OUTDIR/${SAMPLE}_mapped.bam -o $OUTPUT_mapped/${SAMPLE}_mapped_sorted.bam

echo "$SAMPLE extract fastq unmapped"
samtools fastq -@ $NTHREADS $OUTDIR/${SAMPLE}_unmapped_sorted.bam -1 $OUTPUT_unmapped/${SAMPLE}_host_removed_R1.fastq.gz -2 $OUTPUT_unmapped/${SAMPLE}_host_removed_R2.fastq.gz -0 /dev/null -s /dev/null -n

cd $OUTPUT_mapped

echo "$SAMPLE filter mapped by quality"
samtools view -@ $NTHREADS -b -q 30 ${SAMPLE}_mapped_sorted.bam > ${SAMPLE}_qf.bam

echo "$SAMPLE sorting qf"
samtools sort -@ $NTHREADS -n ${SAMPLE}_qf.bam -o ${SAMPLE}_qf2.bam

echo "$SAMPLE add mate score tag"
samtools fixmate -@ $NTHREADS -m -O bam ${SAMPLE}_qf2.bam ${SAMPLE}_mstag.bam

echo "$SAMPLE sorting ms"
samtools sort -@ $NTHREADS ${SAMPLE}_mstag.bam -o ${SAMPLE}_mstag2.bam

echo "$SAMPLE remove duplicates"
samtools markdup -@ $NTHREADS -r ${SAMPLE}_mstag2.bam ${SAMPLE}.bam

echo "$SAMPLE index"
samtools index -@ $NTHREADS ${SAMPLE}.bam

echo "$SAMPLE coverage"
samtools depth ${SAMPLE}.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > $OUTDIR/coverage/${SAMPLE}.txt

rm *_qf.bam *_qf2.bam *_mstag.bam *_mstag2.bam

cd $OUTDIR
rm *.bam
cd $DATADIR
done
```

## HAFpipe

Files `genotypes.vcf` and `SP_combined.fasta` are available [here](https://www.dropbox.com/sh/502iy76q433907o/AABvuUZEUHTZu_RmXp2C_Bfua?dl=0). The output of this script is available inside the `data_HAFpipe` folder and can be used to replicate the analyses below.

```Bash
PRJDIR=path_to_project_directory
DATADIR=path_to_BAM_folder

cd $PRJDIR
wget https://bitbucket.org/dkessner/harp/downloads/harp_linux_140925_103521.zip
unzip harp_linux_140925_103521.zip
rm harp_linux_140925_103521.zip
git clone https://github.com/petrov-lab/HAFpipe-line
mkdir haf_parallel

export PATH=$PATH:$PRJDIR/harp_linux_140925_103521/bin
export PATH=$PATH:$PRJDIR/HAFpipe-line

cd $PRJDIR/haf_parallel

chrArray=("ChrS01" "ChrS02" "ChrS03" "ChrS04" "ChrS05" "ChrS06" "ChrS07" "ChrS08" "ChrS09" "ChrS10" "ChrS11" "ChrS12" "ChrS13" "ChrS14" "ChrS15" "ChrS16" "ChrS17" "ChrS18" "ChrS19" "ChrS20")

function processHAF {
file=$(basename $1 ".bam")
chr=$2
PRJDIR=path_to_project_directory
DATADIR=path_to_BAM_folder
REFDIR=path_to_ref_genome
OUTDIR=$PRJDIR/haf_parallel
echo "Processing sample $file chromosome $chr"
HAFpipe_wrapper.sh -t 1,2,3,4 \
-v $REFDIR/genotypes.vcf \
-c $chr \
-s $OUTDIR/${file}.${chr}.snpTable \
-i simpute \
-b $DATADIR/${file}.bam \
-e sanger \
-a 0 \
-w 100 \
-g 30 \
-r $REFDIR/SP_combined.fasta \
-o $OUTDIR
}

export -f processHAF

parallel -j 16 processHAF  ::: $DATADIR/*.bam ::: "${chrArray[@]}"

find . -type f ! -name "*.freqs" -exec rm -rf {} \; 

for file in $DATADIR/*.bam
do
SAMPLE=$(basename $file ".bam")
cat ${SAMPLE}.bam.* > ${SAMPLE}.txt
done 
```

## Kraken2

We used `Kraken2` to characterize the plant-associated microbial community using the pool-seq data that did not match the *S. polyrhiza* reference genome. The output are available inside the folder `kraken_bacteria_fungi` and `kraken_algae` for replicating analyses.

### Bacteria and fungi

```Bash
PRJDIR=path_to_project_director
DATADIR=path_to_Bowtie_output
DBDIR=$PRJDIR/krakenDB
OUTDIR_kraken=$PRJDIR/kraken
OUTDIR_bracken=$PRJDIR/bracken
NTHREADS=32

cd $DBDIR

wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20210517.tar.gz
tar -xvzf k2_pluspf_20210517.tar.gz

cd $DATADIR

for file in *_host_removed_R1.fastq.gz
do
SAMPLE=$(basename $file "_host_removed_R1.fastq.gz")
cd $OUTDIR
kraken2 --use-names --threads $NTHREADS \
--db $DBDIR \
--gzip-compressed \
--paired $DATADIR/${SAMPLE}_host_removed_R1.fastq.gz \
$DATADIR/${SAMPLE}_host_removed_R2.fastq.gz \
--output ${SAMPLE}.txt --report ${SAMPLE}.report.txt
done 

cd $OUTDIR_kraken

for file in *.report.txt
do
SAMPLE=$(basename $file ".report.txt")
bracken -d $DBDIR -i ${SAMPLE}.report.txt -o $OUTDIR_bracken/${SAMPLE}.txt -l S
done 
```

Using R, we process the data and produce a `phyloseq` object.

```R
library("phyloseq")

files <- dir(path = "data", pattern = "*.txt")
data <- sapply(files, function(x){
  raw.data <- read.table(paste0("data/", x), sep="\t", header = T)[,c("name", "new_est_reads")]
  colnames(raw.data)[2] = x
  return(raw.data)},
  simplify = FALSE,USE.NAMES = TRUE)
data <- Reduce(function(x,y) {merge(x,y,all=T, by="name")}, data)
data[is.na(data)]<-0
colnames(data) <- sub('\\.[^.]+$', '', colnames(data))

row.names(data) <- paste('Taxon', 1:nrow(data), sep = '_')
taxa <- as.data.frame(data[,c(1)])
row.names(taxa) <- row.names(data) 
colnames(taxa)[1] <- "Genus"
write.table(taxa, "taxa_bacteria_fungi.txt")

metadata <- read.table("metadata.txt", header = T, sep = "\t", row.names = 1)
metadata <- sample_data(metadata)

otutable <- as.matrix(data[,c(-1)])
otutable <- otu_table(otutable, taxa_are_rows=TRUE)

colnames(otutable) %in% rownames(metadata)

GM <- merge_phyloseq(otutable, metadata)

GM <- filter_taxa(GM, function (x) {sum(x > 0) > 1}, prune=TRUE)
GM


GM.kracken <- GM
gdata::keep(GM.kracken, sure = TRUE)
save.image(file = 'kraken_bacteria_fungi.rds')
```

For replicating analyses, the files `kraken_bacteria_fungi.rds` and `taxa_bacteria_fungi.txt` are available inside the `data_kraken` folder, while `metadata.txt` is available inside the `reference` folder.

### Algae

First, we create a reference dataset. The file `Acc_List.txt` with the accession IDs of the selected reference genomes is available inside the `reference` folder.

```Bash
PRJDIR=path_to_project_director
DATADIR=path_to_Bowtie_output
DBDIR=$PRJDIR/krakenDB
OUTDIR=$PRJDIR/kraken
OUTDIR_bracken=$PRJDIR/bracken
NTHREADS=32

mkdir ref_database_algae
cd ref_database_algae
mkdir fasta

cat Acc_List.txt | while read line || [[ -n $line ]];
do
echo "Downloading entry $line"
esearch -db nucleotide -query $line < /dev/null | efetch -format fasta > ./fasta/${line}.fasta
done

kraken2-build --download-taxonomy --db krakenDBcustom
kraken2-build --download-library plant --db krakenDBcustom

for file in ./fasta/*.fasta
do
kraken2-build --add-to-library $file --db krakenDBcustom
done

kraken2-build --build --db krakenDBcustom --threads 32

bracken-build -d krakenDBcustom -t 32
```

Then, we re-run Kraken2.

```Bash
PRJDIR=path_to_project_director
DATADIR=path_to_Bowtie_output
DBDIR=$PRJDIR/ref_database_algae
OUTDIR_kraken=$PRJDIR/kraken
OUTDIR_bracken=$PRJDIR/bracken
NTHREADS=32

cd $DATADIR

for file in *_host_removed_R1.fastq.gz
do
SAMPLE=$(basename $file "_host_removed_R1.fastq.gz")
cd $OUTDIR
kraken2 --use-names --threads $NTHREADS \
--db $DBDIR \
--gzip-compressed \
--paired $DATADIR/${SAMPLE}_host_removed_R1.fastq.gz \
$DATADIR/${SAMPLE}_host_removed_R2.fastq.gz \
--output ${SAMPLE}.txt --report ${SAMPLE}.report.txt
done 

cd $OUTDIR_kraken

for file in *.report.txt
do
SAMPLE=$(basename $file ".report.txt")
bracken -d $DBDIR -i ${SAMPLE}.report.txt -o $OUTDIR_bracken/${SAMPLE}.txt -l S
done 
```

Using R, we process the data and produce a `phyloseq` object.

```R
library("phyloseq")

files <- dir(path = "data", pattern = "*.txt")
data <- sapply(files, function(x){
  raw.data <- read.table(paste0("data/", x), sep="\t", header = T)[,c("name", "new_est_reads")]
  colnames(raw.data)[2] = x
  return(raw.data)},
  simplify = FALSE,USE.NAMES = TRUE)
data <- Reduce(function(x,y) {merge(x,y,all=T, by="name")}, data)
data[is.na(data)]<-0
colnames(data) <- sub('\\.[^.]+$', '', colnames(data))

row.names(data) <- paste('Taxon', 1:nrow(data), sep = '_')
taxa <- as.data.frame(data[,c(1)])
row.names(taxa) <- row.names(data) 
colnames(taxa)[1] <- "Genus"
write.table(taxa, "taxa_algae.txt")

metadata <- read.table("metadata.txt", header = T, sep = "\t", row.names = 1)
metadata <- sample_data(metadata)

otutable <- as.matrix(data[,c(-1)])
otutable <- otu_table(otutable, taxa_are_rows=TRUE)

colnames(otutable) %in% rownames(metadata)

GM <- merge_phyloseq(otutable, metadata)

GM <- filter_taxa(GM, function (x) {sum(x > 0) > 1}, prune=TRUE)
GM


GM.kracken <- GM
gdata::keep(GM.kracken, sure = TRUE)
save.image(file = 'kraken_algae.rds')
```

For replicating analyses, the files `kraken_algae.rds` and `taxa_algae.txt` are available inside the `data_kraken` folder, while `metadata.txt` is available inside the `reference` folder.

## Figure 3

Load libraries

```R
library("dplyr")
library("reshape")
library("ggplot2")
library("lme4")
library("car")
library("emmeans")
library("cowplot")
library("ggpubr")
```

Load data and estimate genotype frequency per sample.

```R
single.gen.fun <- function(x, y, z){
  df1 <- read.table(x)
  df2 <- read.table(y)
  df3 <- read.table(z)
  colnames(df1) <- c("position", "begin", "end", "Sp21", "Sp56", "Sp58", "Sp65")
  colnames(df2) <- c("position", "begin", "end", "Sp21", "Sp56", "Sp58", "Sp65")
  colnames(df3) <- c("position", "begin", "end", "Sp21", "Sp56", "Sp58", "Sp65")
  df.genotype <- df1
  df.genotype$Sp21 <- rowMeans(cbind(df1$Sp21, df2$Sp21, df3$Sp21))
  df.genotype$Sp56 <- rowMeans(cbind(df1$Sp56, df2$Sp56, df3$Sp56))
  df.genotype$Sp58 <- rowMeans(cbind(df1$Sp58, df2$Sp58, df3$Sp58))
  df.genotype$Sp65 <- rowMeans(cbind(df1$Sp65, df2$Sp65, df3$Sp65))
  df.genotype <- df.genotype %>% filter_at(vars(4:7), any_vars(. > 0.90))
  df.genotype$position <- paste(df.genotype$position, df.genotype$begin, df.genotype$end, sep = "_")
  df.genotype$begin <- NULL
  df.genotype$end <- NULL
  return(df.genotype)
}

df.Sp21 <- single.gen.fun("data_HAFpipe/A_7.txt", "data_HAFpipe/A_8.txt", "data_HAFpipe/A_9.txt")
df.Sp56 <- single.gen.fun("data_HAFpipe/A_10.txt", "data_HAFpipe/A_11.txt", "data_HAFpipe/A_12.txt")
df.Sp58 <- single.gen.fun("data_HAFpipe/A_13.txt", "data_HAFpipe/A_14.txt", "data_HAFpipe/A_15.txt")
df.Sp65 <- single.gen.fun("data_HAFpipe/A_16.txt", "data_HAFpipe/A_17.txt", "data_HAFpipe/A_18.txt")

comvals <- Reduce(intersect, list(df.Sp21$position, df.Sp56$position, df.Sp58$position, df.Sp65$position))

files <- dir(path = "data_HAFpipe", pattern = "*.txt")

data <- sapply(files, function(x){
  raw.data <- read.table(paste0("data_HAFpipe/", x))
  colnames(raw.data) <- c("position", "begin", "end", "Sp21", "Sp56", "Sp58", "Sp65")
  raw.data$position <- paste(raw.data$position, raw.data$begin, raw.data$end, sep = "_")
  raw.data$sample <- paste0(x)
  raw.data$sample  <- sub('\\.[^.]+$', '', raw.data$sample)
  raw.data$begin <- NULL
  raw.data$end <- NULL
  return(raw.data)},
  simplify = FALSE,USE.NAMES = TRUE)

data <- Reduce(function(x,y) {rbind(x,y)}, data)
data <- data[which(data$position %in% comvals),]
data$position <- NULL

df <- data %>% group_by(sample) %>% summarize(Sp21 = mean(Sp21),
                                              Sp56 = mean(Sp56),
                                              Sp58 = mean(Sp58),
                                              Sp65 = mean(Sp65))

metadata <- read.table("reference/metadata.txt", header = T)

df2 <- merge(metadata, df, by.x = "SampleID", by.y = "sample")
colnames(df2)
mdata <- melt(df2, id=c("SampleID", "category", "pond", "cage", "week")) 
```

Correct genotype frequency using data from the artificial mix

```R
mdata4 <- mdata[which(mdata$category == "mix"),]
mdata4 <- mdata4 %>% group_by(cage, week, variable) %>% summarise(value = mean(value))
mdata4 <- mdata4[which(mdata4$cage == "mix_pond"),]
mdata4 <- mdata4 %>% mutate(value = 0.25/value)

mdata0 <- mdata[which(mdata$cage == "mix_pond"),]
mdata0 <- merge(mdata0, mdata4, by = "variable") %>% mutate(value = value.x * value.y)
mdata0 <- mdata0 %>% select(SampleID, pond, cage.x, week.x, variable, value)
colnames(mdata0) <- c("sampleID", "pond", "cage", "week", "genotype", "value")

mdata2 <- mdata[which(mdata$category == "pond"),]
mdata2 <- merge(mdata2, mdata4, by = "variable") %>% mutate(value = value.x * value.y)
mdata2 <- mdata2 %>% select(SampleID, pond, cage.x, week.x, variable, value)
colnames(mdata2) <- c("sampleID", "pond", "cage", "week", "genotype", "value")

mdata.df <- rbind(mdata0, mdata2)
```

Plot the data

```R
plot.fun <- function(gtype, ylabel, siglabel){
  mdata.plot <- mdata.df %>% filter(genotype == gtype)
  px <- ggplot(mdata.plot, aes(x = week, y = value, fill = cage, color = cage)) +
      theme_bw(base_size = 14) +
      stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + sd(x), fun.min = function(x) mean(x) - sd(x), position=position_dodge(0)) +
      theme(axis.title.x=element_blank(),
            axis.text.x = element_text(color="black"),
            axis.text.y = element_text(color="black"),
            axis.title.y = element_text(color="black"),
            legend.position = "none",
            panel.grid = element_blank()) +
      scale_fill_manual(name = "Legend", values=c("#252525", "#4daf4a", "#e41a1c"), labels = c("Starting", "Control", "Herbivory"), breaks = c("mix_pond", "A_ctr", "B_herb")) +
      scale_color_manual(name = "Legend", values=c("#252525", "#4daf4a", "#e41a1c"), labels = c("Starting", "Control", "Herbivory"), breaks = c("mix_pond", "A_ctr", "B_herb")) +
      scale_x_discrete(name = "Legend", labels = c("0", "8", "12"), breaks = c("w0", "w08", "w12")) +
    
      geom_segment(aes(x=1,xend=2,
                       y= mean(mdata.plot[which(mdata.plot$cage == "mix_pond" & mdata.plot$week == "w0"),]$value),
                       yend=mean(mdata.plot[which(mdata.plot$cage == "A_ctr" & mdata.plot$week == "w08"),]$value)),
                   color = "#4daf4a") +
      
      geom_segment(aes(x=2,xend=3,
                       y= mean(mdata.plot[which(mdata.plot$cage == "A_ctr" & mdata.plot$week == "w08"),]$value),
                       yend=mean(mdata.plot[which(mdata.plot$cage == "A_ctr" & mdata.plot$week == "w12"),]$value)),
                   color = "#4daf4a") +
      
      geom_segment(aes(x=1,xend=2,
                       y= mean(mdata.plot[which(mdata.plot$cage == "mix_pond" & mdata.plot$week == "w0"),]$value),
                       yend=mean(mdata.plot[which(mdata.plot$cage == "B_herb" & mdata.plot$week == "w08"),]$value)),
                   color = "#e41a1c") +
      
      geom_segment(aes(x=2,xend=3,
                       y= mean(mdata.plot[which(mdata.plot$cage == "B_herb" & mdata.plot$week == "w08"),]$value),
                       yend=mean(mdata.plot[which(mdata.plot$cage == "B_herb" & mdata.plot$week == "w12"),]$value)),
                   color = "#e41a1c") +
      
      annotate("text", x = 3.3, 
               y =  (mean(mdata.plot[which(mdata.plot$cage == "A_ctr" & mdata.plot$week == "w12"),]$value) +
                              mean(mdata.plot[which(mdata.plot$cage == "B_herb" & mdata.plot$week == "w12"),]$value))/2,
               label = paste0(siglabel), size = 5) +
      
      ggtitle(paste0(gtype)) +
      
      ylim(0,0.6) +
      labs(y = paste0(ylabel))
  
  return(px)
}

px1 <- plot.fun("Sp21", "Genotype frequency", " ")
px2 <- plot.fun("Sp56", " ", " ")
px3 <- plot.fun("Sp58", "Genotype frequency", " ")
px4 <- plot.fun("Sp65", " ", " ")


px <- ggarrange(px1, px2, px3, px4,
                ncol = 2, nrow = 2,  labels = c("a.", "b.", "c.", "d."), common.legend = T)

px
```

## Table S3

```R
mdata.df2 <- mdata.df[which(mdata.df$cage != "mix_pond"),]
model <- lmer(value ~ cage * genotype * (1|week) * (1|pond), data = mdata.df2)
Anova(model)
m1 <- emmeans(model, "cage" , by = c("genotype"))
pairs(m1)
```

## Table S5

Load libraries

```R
library("purrr")
library("ape")
library("car")
library("data.table")
library("DESeq2")
library("dplyr")
library("emmeans")
library("ggplot2")
library("ggpmisc")
library("ggpubr")
library("grid")
library("lme4")
library("phyloseq")
library("picante")
library("plyr")
library("Rmisc")
library("tidyr")
library("tibble")
library("vegan")
library("Hmisc")
library("limma")
library("scales")
library("RVAideMemoire")
```

Bacteria and fungi. Load data, normalize, and run PERMANOVA.

```R
load(file = 'data_kraken/kraken_bacteria_fungi.rds')
kraken.tax <- read.table("data_kraken/taxa_bacteria_fungi.txt", header = T, row.names = NULL)

GM.kracken <- subset_samples(GM.kracken, category =="pond")
GM.kracken <- filter_taxa(GM.kracken, function (x) {sum(x > 0) > 1}, prune=TRUE)

norm.deseq <- function(x){
  diagdds = phyloseq_to_deseq2(x, ~ 1)
  ts = counts(diagdds)
  geoMeans = apply(ts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  diagdds = estimateSizeFactors(diagdds, geoMeans=geoMeans)
  diagdds = estimateDispersions(diagdds)
  diagvst = getVarianceStabilizedData(diagdds)
  diagdds.c <- removeBatchEffect(diagvst)
  diagdds.c[diagdds.c<0] <- 0
  ps.16Sn <- x
  otu_table(ps.16Sn) <- otu_table(diagdds.c, taxa_are_rows = TRUE)
}

GM <- norm.deseq(GM.kracken)
sampledf <- data.frame(sample_data(GM.kracken))
dist.mat <- phyloseq::distance(GM, method = "bray")
perm <- how(nperm = 999)
set.seed(100)
setBlocks(perm) <- with(sampledf, pond)
pmv <- adonis2(dist.mat ~  cage * week, data = sampledf, permutations = perm)
pmv
```

Algae. Load data, normalize, and run PERMANOVA.

```R
load(file = 'data_kraken/kraken_algae.rds')
kraken.tax <- read.table("data_kraken/taxa_algae.txt", header = T, row.names = NULL)

GM.kracken <- subset_samples(GM.kracken, category =="pond")
GM.kracken <- filter_taxa(GM.kracken, function (x) {sum(x > 0) > 1}, prune=TRUE)

norm.deseq <- function(x){
  diagdds = phyloseq_to_deseq2(x, ~ 1)
  ts = counts(diagdds)
  geoMeans = apply(ts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  diagdds = estimateSizeFactors(diagdds, geoMeans=geoMeans)
  diagdds = estimateDispersions(diagdds)
  diagvst = getVarianceStabilizedData(diagdds)
  diagdds.c <- removeBatchEffect(diagvst)
  diagdds.c[diagdds.c<0] <- 0
  ps.16Sn <- x
  otu_table(ps.16Sn) <- otu_table(diagdds.c, taxa_are_rows = TRUE)
}

GM <- norm.deseq(GM.kracken)
sampledf <- data.frame(sample_data(GM.kracken))
dist.mat <- phyloseq::distance(GM, method = "bray")
perm <- how(nperm = 999)
set.seed(100)
setBlocks(perm) <- with(sampledf, pond)
pmv <- adonis2(dist.mat ~  cage * week, data = sampledf, permutations = perm)
pmv
```
