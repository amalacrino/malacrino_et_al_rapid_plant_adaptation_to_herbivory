# Indoors assays

All the data to replicate analyses below is available inside the `data` folder.

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
## Figure S4

```R
df1 <- read.table("data/GA_indoors.txt", sep = "\t", header = T)
df2 <- read.table("data/HA_indoors.txt", sep = "\t", header = T)

df1$rel.count <- (log(df1$counts) - log(30))/13

px11 <- ggplot(df1, aes(x = genotype, y = rel.count, fill = genotype)) +
  theme_bw(base_size = 14) +
  geom_violin(trim=FALSE) +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + sd(x), fun.min = function(x) mean(x) - sd(x), position=position_dodge(0.9), color="black") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        legend.position = "none") +
  scale_fill_manual(name = "Legend", values=c("#a6611a", "#dfc27d", "#80cdc1", "#018571"), labels = c("Sp21", "Sp56", "Sp58", "Sp65"), breaks = c("Sp21", "Sp56", "Sp58", "Sp65")) +
  ylim(0,0.3) +
  labs(y = "Relative growth rate")

px12 <- ggplot(df2, aes(x = genotype, y = consumed, fill = genotype)) +
  theme_bw(base_size = 14) +
  geom_violin(trim=FALSE) +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + sd(x), fun.min = function(x) mean(x) - sd(x), position=position_dodge(0.9), color="black") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        legend.position = "none") +
  scale_fill_manual(name = "Legend", values=c("#a6611a", "#dfc27d", "#80cdc1", "#018571"), labels = c("Sp21", "Sp56", "Sp58", "Sp65"), breaks = c("Sp21", "Sp56", "Sp58", "Sp65")) +
  ylim(0,70) +
  labs(y = "# of consumed fronds")

px <- ggarrange(px11, px12, ncol = 2, labels = c("a.", "b."))
px

model <- lm(rel.count ~ genotype, data = df1)
Anova(model)
m1 <- emmeans(model, "genotype")
pairs(m1)

model <- lm(consumed ~ genotype, data = df2)
Anova(model)
m1 <- emmeans(model, "genotype")
pairs(m1)
```

## Figure S5

```R
df <- read.table("data/GA_indoors_SG.txt", sep = "\t", header = T)

df$countChange <- (log(df$count_t13) - log(200))/13

px8 <- ggplot(df, aes(x = genotype, y = countChange, fill = treatment)) +
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
  labs(y = "Relative growth rate")

px9 <- ggplot(df, aes(x = genotype, y = area_t13, fill = treatment)) +
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
  labs(y = "Area per frond (cm2)")

px10 <- ggplot(df, aes(x = genotype, y = biomass, fill = treatment)) +
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
  labs(y = "Weight per frond (mg)")

px <- ggarrange(px8, px9, px10, ncol = 3, labels = c("a.", "b.", "c."))
px
```

## Table S4

```R
model <- lmer(countChange ~ genotype * treatment * (1|pond), data = df)
Anova(model)

model <- lmer(area_t13 ~ genotype * treatment * (1|pond), data = df)
Anova(model)

model <- lmer(biomass ~ genotype * treatment * (1|pond), data = df)
Anova(model)
```

## Figure 4

```R
df <- read.table("data/SynPop.txt", sep = "\t", header = T)

px5 <- ggplot(df, aes(x = genotype, y = consumed, fill = treatment)) +
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
  ylim(-40, 130) +
  labs(y = "# of consumed fronds")

px5

model <- lmer(consumed ~ genotype * treatment * (1|pond), data = df)
Anova(model)
m1 <- emmeans(model, "treatment", by = "genotype")
pairs(m1)
```
