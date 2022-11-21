# Induced responses contribute to rapid plant adaptation to herbivory outdoors 

**Antonino Malacrinò, Laura Böttner, Sara Nouere, Meret Huber, Martin Schäfer, Shuqing Xu**

## Abstract
Herbivory-induced responses in plants are a typical example of phenotypic plasticity and are thought to be important for adaptation to herbivory. However, direct evidence demonstrating that induced responses increase the multigenerational fitness of plants remains scarce. Here, we experimentally evolved populations of an aquatic plant and its native herbivore for more than 30 generations outdoors. We found that herbivory rapidly increased plant resistance via altering the genotype frequencies, herbivory-induced phenotypic plasticity, and plant microbiota. Together, this study provides direct evidence for the adaptive role of herbivory-induced responses in plants.

# Disclaimer

This repository contains the main components used to process the raw data and to analyze it. Raw data is available at NCBI SRA under the BioProject numbers `PRJNA893536` (16S amplicon metagenomics), `PRJNA893537` (18S amplicon metagenomics) and `PRJNA849359` (pool-seq).

Our pipeline included:
* TrimGalore (Krueger et al. [2021](https://doi.org/10.5281/zenodo.5127899))
* Bowtie2 (Langmead et al. [2019](https://academic.oup.com/bioinformatics/article/35/3/421/5055585))
* Samtools (Danecek et al. [2021](https://doi.org/10.1093/gigascience/giab008))
* HAFpipe (Tilk et al. [2019](https://academic.oup.com/g3journal/article/9/12/4159/6028100)) and HARP (Kessner et al. [2013](https://academic.oup.com/mbe/article/30/5/1145/993804))
* Kraken2 (Wood et al. [2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0))
* DADA2 (Callahan et al. [2016](https://www.nature.com/articles/nmeth.3869))
* R  (R Core Team [2022](https://www.R-project.org/))

# Code

### **1.** [Duckweed monitoring outdoors](/1_monitoring)
Data and code for reproducing Fig. 1 and Tab. S2.

### **2.** [Growth/Herbivory assays outdoors](/2_GH_assays)
Data and code for reproducing Fig. 2, Fig. 5a, Fig. S2, Tab. S6 and Tab. S7.

### **3.** [Pool-seq](/3_pool_seq)
Data and code for reproducing Fig. 3, Tab. S3, and Tab. S5.

### **4.** [Indoors assays](/4_indoors)
Data and code for reproducing Fig. 4, Fig. S4, Fig. S5, and Tab. S4.

### **5.** [Amplicon metagenomics](/5_amplicon_metagenomics)
Data and code for reproducing Fig. 5b, Fig. 5c, and Tab. S8.

### **6.** [Metabolomics](/6_metabolomics)
Data and code for reproducing Fig. S3.

### **7.** [Water nutrients](/7_nutrients)
Data and code for reproducing Fig. S6 and Tab. S9.
