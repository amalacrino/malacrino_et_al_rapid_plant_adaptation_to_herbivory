# Metabolomics

### Load packages

```R
library("dplyr")
library("lme4")
library("car")
library("emmeans")
library("ggplot2")
library("cowplot")
library("vegan")
library("tidyr")
library("data.table")
library("ggrepel")
```

### Process 2021 data

```R
data <- read.table("data/data_metabolome_2021.txt", sep = "\t", header = T)

data.long <- data %>% gather(metabolite, value, 7:55) 
group.data.long  <- group_by(data.long , metabolite)
list.mol <-unique(c(as.character(data.long$metabolite)))

model_calculator.1 <- sapply(list.mol,  
                             function(x){
                               data.s <- group.data.long[which(group.data.long$metabolite==x),]
                               model <- lmer(value ~ cage * (1|pond) * (1|group), data = data.s)
                               aaa <-  Anova(model)
                               aaa$sig = c(rep('',length(aaa$`Pr(>Chisq)`)))
                               makeStars <- function(x){
                                 stars <- c("****", "***", "**", "*", "ns")
                                 vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
                                 i <- findInterval(x, vec)
                                 stars[i] }
                               aaa$sig <- makeStars(aaa$`Pr(>Chisq)`)
                               aaa <- aaa[1,]
                               return(aaa)},
                             simplify = FALSE,USE.NAMES = TRUE)

res <- do.call(rbind, model_calculator.1)
res <- setDT(res, keep.rownames = TRUE)[]

table.analysis.1 <- data.long %>% 
                      filter(cage == "control") %>%
                      group_by(metabolite) %>% 
                      dplyr::summarize(control = mean(value, na.rm=TRUE))

table.analysis.2 <- data.long %>% 
                      filter(cage == "herbivory") %>%
                      group_by(metabolite) %>% 
                      dplyr::summarize(herbivory = mean(value, na.rm=TRUE))


res.def <- merge(res, table.analysis.1, by.x = "rn", by.y = "metabolite")
res.def <- merge(res.def, table.analysis.2, by.x = "rn", by.y = "metabolite")

res.def$log2FoldChange <-  log2(res.def$herbivory/res.def$control)

de <- res.def
names(de)[names(de) == 'Pr(>Chisq)'] <- 'pvalue'
de$diffexpressed <- "no changes"
de$diffexpressed[de$log2FoldChange > 0.1 & de$pvalue < 0.05] <- "herbivory"
de$diffexpressed[de$log2FoldChange < -0.1 & de$pvalue < 0.05] <- "control"

plot <- ggplot(data=de) +
        theme_bw(base_size = 12) +
        geom_point(aes(x = log2FoldChange, y = -log10(pvalue), colour = diffexpressed)) +
        geom_text_repel(aes(x = log2FoldChange, y = -log10(pvalue), label = ifelse(diffexpressed != "no changes", rn,"")), max.overlaps = 20) +
        scale_color_manual(values=c("blue", "red", "black")) +
        geom_vline(xintercept=0, col="black", linetype = "longdash") +
        geom_hline(yintercept=-log10(0.05), col="black", linetype = "longdash") +
        theme(legend.justification=c(0.01,0.99), legend.position=c(0.01,0.99),
              panel.grid = element_blank(),
              legend.title=element_blank(), 
              legend.background = element_rect(color = NA),
              legend.key = element_rect(color = NA),
              panel.border = element_rect(colour = "black", fill=NA, size=1)) +
        xlab(expression(paste(Log[2], " Fold Changes"))) +
        ylab(expression(paste(-Log[10], " P"))) +
        scale_color_manual(name = "Legend", values=c("#4daf4a", "#e41a1c", "#000000"), labels = c("Control", "Herbivory", "No changes"), breaks = c("control", "herbivory", "no changes"))
        
plot
```

### Process 2022 data

```R
data <- read.table("data/data_metabolome_2022.txt", sep = "\t", header = T)

data.long <- data %>% gather(metabolite, value, 5:55) 
group.data.long  <- group_by(data.long , metabolite)
list.mol <-unique(c(as.character(data.long$metabolite)))

model_calculator.1 <- sapply(list.mol,  
                             function(x){
                               data.s <- group.data.long[which(group.data.long$metabolite==x),]
                               model <- lmer(value ~ cage * (1|pond), data = data.s)
                               aaa <-  Anova(model)
                               aaa$sig = c(rep('',length(aaa$`Pr(>Chisq)`)))
                               makeStars <- function(x){
                                 stars <- c("****", "***", "**", "*", "ns")
                                 vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
                                 i <- findInterval(x, vec)
                                 stars[i] }
                               aaa$sig <- makeStars(aaa$`Pr(>Chisq)`)
                               aaa <- aaa[1,]
                               return(aaa)},
                             simplify = FALSE,USE.NAMES = TRUE)

res <- do.call(rbind, model_calculator.1)
res <- setDT(res, keep.rownames = TRUE)[]

table.analysis.1 <- data.long %>% 
                      filter(cage == "control") %>%
                      group_by(metabolite) %>% 
                      dplyr::summarize(control = mean(value, na.rm=TRUE))

table.analysis.2 <- data.long %>% 
                      filter(cage == "herbivory") %>%
                      group_by(metabolite) %>% 
                      dplyr::summarize(herbivory = mean(value, na.rm=TRUE))


res.def <- merge(res, table.analysis.1, by.x = "rn", by.y = "metabolite")
res.def <- merge(res.def, table.analysis.2, by.x = "rn", by.y = "metabolite")

res.def$log2FoldChange <-  log2(res.def$herbivory/res.def$control)

de <- res.def
names(de)[names(de) == 'Pr(>Chisq)'] <- 'pvalue'
de$diffexpressed <- "no changes"
de$diffexpressed[de$log2FoldChange > 0.1 & de$pvalue < 0.05] <- "herbivory"
de$diffexpressed[de$log2FoldChange < -0.1 & de$pvalue < 0.05] <- "control"

plot <- ggplot(data=de) +
        theme_bw(base_size = 12) +
        geom_point(aes(x = log2FoldChange, y = -log10(pvalue), colour = diffexpressed)) +
        geom_text_repel(aes(x = log2FoldChange, y = -log10(pvalue), label = ifelse(diffexpressed != "no changes", rn,"")), max.overlaps = 20) +
        scale_color_manual(values=c("blue", "red", "black")) +
        geom_vline(xintercept=0, col="black", linetype = "longdash") +
        geom_hline(yintercept=-log10(0.05), col="black", linetype = "longdash") +
        theme(legend.justification=c(0.01,0.99), legend.position=c(0.01,0.99),
              panel.grid = element_blank(),
              legend.title=element_blank(), 
              legend.background = element_rect(color = NA),
              legend.key = element_rect(color = NA),
              panel.border = element_rect(colour = "black", fill=NA, size=1)) +
        xlab(expression(paste(Log[2], " Fold Changes"))) +
        ylab(expression(paste(-Log[10], " P"))) +
        scale_color_manual(name = "Legend", values=c("#4daf4a", "#e41a1c", "#000000"), labels = c("Control", "Herbivory", "No changes"), breaks = c("control", "herbivory", "no changes"))
        
plot
```
