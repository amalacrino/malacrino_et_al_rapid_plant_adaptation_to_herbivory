
# Growth and herbivory assays outdoors

All the data to replicate the analyses below is available inside the folder `data`.

## Load libraries

```R
library("dplyr")
library("lme4")
library("car")
library("emmeans")
library("ggplot2")
library("cowplot")
library("ggpubr")
library("data.table")

nperm <- 999
```

## Figure 2

```R
df <- read.table("data/HA_outdoors.txt", sep = "\t", header = T)

set.seed(100)
refs.c <- df %>% filter(round == "reference") %>% filter(cage == "Control")  %>% sample_n(., nperm, replace = T)
refs.h <- df %>% filter(round == "reference") %>% filter(cage == "Herbivory")  %>% sample_n(., nperm, replace = T)
expl.c <- df %>% filter(round != "reference") %>% filter(cage == "Control") %>% sample_n(., nperm, replace = T)
expl.h <- df %>% filter(round != "reference") %>% filter(cage == "Herbivory")  %>% sample_n(., nperm, replace = T)
refs <- rbind(refs.c, refs.h)
expl <- rbind(expl.c, expl.h)

norm <- data.frame(
  round = expl$round,
  pond = expl$pond,
  cage = expl$cage,
  biomass = refs$biomass - expl$biomass,
  area = refs$area - expl$area
)

px5 <- ggplot(df, aes(x = cage, y = consumed, fill = cage)) +
  theme_bw(base_size = 14) +
  geom_violin(trim=FALSE) +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + sd(x), fun.min = function(x) mean(x) - sd(x), position=position_dodge(0.9), color="black") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        legend.position = "none") +
  scale_fill_manual(name = "Legend", values=c("#4daf4a", "#e41a1c"), labels = c("Control", "Herbivory"), breaks = c("Control", "Herbivory")) +
  ylim(-9,65) +
  labs(y = "# of consumed fronds")

px6 <- ggplot(norm, aes(x = cage, y = biomass, fill = cage)) +
  theme_bw(base_size = 14) +
  geom_violin(trim=FALSE) +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + sd(x), fun.min = function(x) mean(x) - sd(x), position=position_dodge(0.9), color="black") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        legend.position = "none") +
  scale_fill_manual(name = "Legend", values=c("#4daf4a", "#e41a1c"), labels = c("Control", "Herbivory"), breaks = c("Control", "Herbivory")) +
  ylim(-10,20) +
  labs(y = "Consumed biomass (mg)")

px7 <- ggplot(norm, aes(x = cage, y = area, fill = cage)) +
  theme_bw(base_size = 14) +
  geom_violin(trim=FALSE) +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + sd(x), fun.min = function(x) mean(x) - sd(x), position=position_dodge(0.9), color="black") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        legend.position = "none") +
  scale_fill_manual(name = "Legend", values=c("#4daf4a", "#e41a1c"), labels = c("Control", "Herbivory"), breaks = c("Control", "Herbivory")) +
  ylim(-10,20) +
  labs(y = "Consumed area (cm2)")

px <- ggarrange(px5, px6, px7,
                ncol = 3, labels = c("a.", "b.", "c."))

px
```

```R
model <- lmer(consumed ~ cage * (1|round) * (1|pond), data = df)
Anova(model)

model <- lmer(biomass ~ cage * (1|round) * (1|pond), data = norm)
Anova(model)

model <- lmer(area ~ cage * (1|round) * (1|pond), data = norm)
Anova(model)
```

## Figure 5a

```R
df <- read.table("data/HA_SG_outdoors.txt", sep = "\t", header = T)

px_counts <- ggplot(df, aes(x = genotype, y = consumed, fill = cage)) +
  theme_bw(base_size = 14) +
  geom_violin(trim=FALSE) +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + sd(x), fun.min = function(x) mean(x) - sd(x), position=position_dodge(0.9), color="black") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        legend.position = "none") +
  scale_fill_manual(name = "Legend", values=c("#4daf4a", "#e41a1c"), labels = c("Control", "Herbivory"), breaks = c("Control", "Herbivory")) +
  ylim(-50,120) +
  labs(y = "# of consumed fronds")
px_counts
```

## Table S7

```R
model <- lmer(consumed ~  cage * (1|pond), data = df[which(df$genotype == "Sp21"),])
Anova(model)

model <- lmer(consumed ~  cage * (1|pond), data = df[which(df$genotype == "Sp56"),])
Anova(model)

model <- lmer(consumed ~  cage * (1|pond), data = df[which(df$genotype == "Sp58"),])
Anova(model)

model <- lmer(consumed ~  cage * (1|pond), data = df[which(df$genotype == "Sp65"),])
Anova(model)
```

## Figure S2

```R
df <- read.table("data/GA_outdoors.txt", sep = "\t", header = T)
df$rel.count <- (log(df$counts) - log(50))/10

px8 <- ggplot(df, aes(x = cage, y = rel.count, fill = cage)) +
  theme_bw(base_size = 14) +
  geom_violin(trim=FALSE) +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + sd(x), fun.min = function(x) mean(x) - sd(x), position=position_dodge(0.9), color="black") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        legend.position = "none") +
  scale_fill_manual(name = "Legend", values=c("#4daf4a", "#e41a1c"), labels = c("Control", "Herbivory"), breaks = c("Control", "Herbivory")) +
  ylim(0,0.2) +
  labs(y = "Rel. growth rate")

refs <- read.table("data/reference_GA_outdoors.txt", sep = "\t", header = T)

set.seed(100)
refs.c <- refs %>% filter(cage == "Control")  %>% sample_n(., nperm, replace = T)
refs.h <- refs %>% filter(cage == "Herbivory")  %>% sample_n(., nperm, replace = T)
expl.c <- df %>% filter(cage == "Control") %>% sample_n(., nperm, replace = T)
expl.h <- df %>% filter(cage == "Herbivory")  %>% sample_n(., nperm, replace = T)
refs <- rbind(refs.c, refs.h)
expl <- rbind(expl.c, expl.h)

norm <- data.frame(
  pond = expl$pond,
  cage = expl$cage,
  biomass = (log(expl$biomass) - log(refs$biomass))/10,
  area = (log(expl$area) - log(refs$area))/10
)

px9 <- ggplot(norm, aes(x = cage, y = biomass, fill = cage)) +
  theme_bw(base_size = 14) +
  geom_violin(trim=FALSE) +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + sd(x), fun.min = function(x) mean(x) - sd(x), position=position_dodge(0.9), color="black") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        legend.position = "none") +
  scale_fill_manual(name = "Legend", values=c("#4daf4a", "#e41a1c"), labels = c("Control", "Herbivory"), breaks = c("Control", "Herbivory")) +
  ylim(0,0.2) +
  labs(y = "Rel. change in biomass")

px10 <- ggplot(norm, aes(x = cage, y = area, fill = cage)) +
  theme_bw(base_size = 14) +
  geom_violin(trim=FALSE) +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + sd(x), fun.min = function(x) mean(x) - sd(x), position=position_dodge(0.9), color="black") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        legend.position = "none") +
  scale_fill_manual(name = "Legend", values=c("#4daf4a", "#e41a1c"), labels = c("Control", "Herbivory"), breaks = c("Control", "Herbivory")) +
  ylim(0,0.8) +
  labs(y = "Rel . change in area")

px <- ggarrange(px8, px9, px10,
                ncol = 3, nrow = 1, labels = c("a.", "b.", "c."))
px

model <- lmer(rel.count ~ cage * (1|pond), data = df)
Anova(model)

model <- lmer(biomass ~ cage * (1|pond), data = norm)
Anova(model)

model <- lmer(area ~ cage * (1|pond), data = norm)
Anova(model)
```

## Table S6

```R
df <- read.table("data/GA_SG_outdoors.txt", sep = "\t", header = T)

list.gen <- (unique(df$genotype))

rgr.cal2 <- function(x){
  df <- df %>% filter(genotype == x)
  w2.c <- df %>% filter(week == "w2") %>% filter(cage == "Control") 
  w2.h <- df %>% filter(week == "w2") %>% filter(cage == "Herbivory") 
  w8.c <- df %>% filter(week == "w8") %>% filter(cage == "Control") 
  w8.h <- df %>% filter(week == "w8") %>% filter(cage == "Herbivory") 
  
  norm1 <- data.frame(genotype = x, treatment = "Control", pond = w8.c$pond, weight = (log(w8.c$biomass) - log(w2.c$biomass))/42, area = (log(w8.c$area) - log(w2.c$area))/42)
  norm2 <- data.frame(genotype = x, treatment = "Herbivory", pond = w8.h$pond, weight = (log(w8.h$biomass) - log(w2.h$biomass))/42, area = (log(w8.h$area) - log(w2.h$area))/42)
  
  norm <- rbind(norm1, norm2)
  return(norm)
}
norm <-lapply(list.gen, rgr.cal2)
norm <- do.call(rbind, norm)

model <- lmer(weight ~ genotype * treatment * (1|pond), data = norm)
Anova(model)

model <- lmer(area ~ genotype * treatment  * (1|pond), data = norm)
Anova(model)

df_ha <- read.table("data/HA_SG_outdoors.txt", sep = "\t", header = T)

norm.cal<- function(x){
  df$area2 <- (df$area/df$count)
  df_ha$area2 <- (df_ha$area/(100-df_ha$consumed))
  
  set.seed(100)
  refs.c <- df %>% filter(genotype == x) %>% filter(cage == "Control") %>% filter(week == "w8") %>% sample_n(., nperm, replace = T)
  refs.h <- df %>% filter(genotype == x) %>% filter(cage == "Herbivory") %>% filter(week == "w8") %>% sample_n(., nperm, replace = T)
  expl.c <- df_ha %>% filter(genotype == x) %>% filter(cage == "Control") %>% sample_n(., nperm, replace = T)
  expl.h <- df_ha %>% filter(genotype == x) %>% filter(cage == "Herbivory")  %>% sample_n(., nperm, replace = T)
  refs <- rbind(refs.c, refs.h)
  expl <- rbind(expl.c, expl.h)
  
  norm <- data.frame(
    genotype = expl$genotype,
    pond = expl$pond,
    cage = expl$cage,
    biomass = refs$biomass - expl$biomass,
    area = refs$area2 - expl$area2
  )
  
  return(norm)
}

norm <-lapply(list.gen, norm.cal)
norm <- do.call(rbind, norm)

model <- lmer(biomass ~ genotype * cage * (1|pond), data = norm)
Anova(model)

model <- lmer(area ~ genotype * cage * (1|pond), data = norm)
Anova(model)
```










