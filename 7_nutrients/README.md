# Water nutrients

### Load packages and data

```R
library("dplyr")
library("ggplot2")
library("cowplot")
library("ggpubr")
library("lme4")
library("car")
library("emmeans")
library("vegan")

df <- read.table("data/data.txt", sep = "\t", header = T)
```

### Table S9

```R
model <- lmer(Cl ~ cage * (1|pond) * (1|week), data = df)
Anova(model)

model <- lmer(NO2 ~ cage * (1|pond) * (1|week), data = df)
Anova(model)

model <- lmer(NO3 ~ cage * (1|pond) * (1|week), data = df)
Anova(model)

model <- lmer(PO4 ~ cage * (1|pond) * (1|week), data = df)
Anova(model)

model <- lmer(SO4 ~ cage * (1|pond) * (1|week), data = df)
Anova(model)

model <- lmer(Na ~ cage * (1|pond) * (1|week), data = df)
Anova(model)

model <- lmer(NH4 ~ cage * (1|pond) * (1|week), data = df)
Anova(model)

model <- lmer(K ~ cage * (1|pond) * (1|week), data = df)
Anova(model)

model <- lmer(Ca ~ cage * (1|pond) * (1|week), data = df)
Anova(model)

model <- lmer(Mg ~ cage * (1|pond) * (1|week), data = df)
Anova(model)

model <- lmer(pH ~ cage * (1|pond) * (1|week), data = df)
Anova(model)
```

### Figure S6 - Plot

```R
df <- df[which(df$week != "week00"),]
df <- df[which(df$week != "week02"),]

df1 <- df[,c(7:16)]
vare.mds <- metaMDS(df1) 
data.scores <- as.data.frame(scores(vare.mds)$sites)
data.scores <- cbind(df[,c(1:6)], data.scores)

nmds <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2, fill = treat)) +
  theme_bw(base_size = 14) +
  stat_ellipse(aes(fill = treat), alpha = 0.2, geom = "polygon", size=0.3, show.legend=T, colour = "black") +
  geom_point(aes(color = treat , fill = treat), size = 5) +
  theme(legend.title=element_blank(), 
        legend.background = element_rect(fill="transparent"),
        legend.key = element_rect(fill="transparent"),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        panel.grid = element_blank()) +
  scale_fill_manual(name = "Legend", values=c("#4daf4a", "#e41a1c"), labels = c("control", "herbivory"), breaks = c("control", "herbivory")) +
  scale_color_manual(name = "Legend", values=c("#4daf4a", "#e41a1c"), labels = c("control", "herbivory"), breaks = c("control", "herbivory"))
nmds
```

### Figure S6 - PERMANOVA

```R
perm <- how(nperm = 1000)
setBlocks(perm) <- with(df, pond, week)
permanova <- adonis2(df1 ~ treat, data = df, permutations = perm, method = "bray")
permanova
```
