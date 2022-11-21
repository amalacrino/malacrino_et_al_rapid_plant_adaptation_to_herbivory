# Duckweed monitoring outdoors.

The dataset for reproducing Figure 1 is available inside the `data` folder.

## Figure 1

```R
library("dplyr")
library("ggplot2")
library("cowplot")
library("ggpubr")
library("lme4")
library("car")
library("emmeans")

df <- read.table("data/data.txt", sep = "\t", header = T)

mdata.plot <- df %>% group_by(year, timepoint, treatment) %>% summarise(area_av = mean(area),
                                                                        area_se = sd(area),
                                                                        weight_av = mean(weight_fr),
                                                                        weight_se = sd(weight_fr))


px1 <- ggplot(mdata.plot, aes(x=timepoint, y=weight_av, color=treatment)) +
  theme_bw(base_size = 14) +
  geom_pointrange(aes(ymin=weight_av-weight_se, ymax=weight_av+weight_se), position = position_dodge(0.9)) +
  theme(axis.title.x = element_text(color="black"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))  +
  scale_color_manual(name = "Legend", values=c("#252525", "#4daf4a", "#e41a1c"), labels = c("starting", "Control", "Herbivory"), breaks = c("Starting", "Control", "Herbivory")) +
  scale_x_continuous(breaks = seq(0, 60, by = 4)) +
  expand_limits(x = c(0,61), y = 0) +

  geom_rect(data=NULL,aes(xmin = 16, xmax = 40 ,ymin=-Inf,ymax=Inf), fill="gray", color = "transparent") +
  
  geom_segment(aes(x=4-0.2, xend=8-0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 4),]$weight_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 8),]$weight_av),
                   color = "#4daf4a") +
  geom_segment(aes(x=8-0.2, xend=12-0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 8),]$weight_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 12),]$weight_av),
               color = "#4daf4a") +

  geom_segment(aes(x=4+0.2, xend=8+0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 4),]$weight_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 8),]$weight_av),
               color = "#e41a1c") +
  geom_segment(aes(x=8+0.2, xend=12+0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 8),]$weight_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 12),]$weight_av),
               color = "#e41a1c") +
  
  geom_segment(aes(x=42-0.2, xend=46-0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 42),]$weight_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 46),]$weight_av),
               color = "#4daf4a") +
  geom_segment(aes(x=46-0.2, xend=50-0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 46),]$weight_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 50),]$weight_av),
               color = "#4daf4a") +
  geom_segment(aes(x=50-0.2, xend=55-0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 50),]$weight_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 55),]$weight_av),
               color = "#4daf4a") +
  geom_segment(aes(x=55-0.2, xend=58-0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 55),]$weight_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 58),]$weight_av),
               color = "#4daf4a") +
  
  geom_segment(aes(x=42+0.2, xend=46+0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 42),]$weight_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 46),]$weight_av),
               color = "#e41a1c") +
  geom_segment(aes(x=46+0.2, xend=50+0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 46),]$weight_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 50),]$weight_av),
               color = "#e41a1c") +
  geom_segment(aes(x=50+0.2, xend=55+0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 50),]$weight_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 55),]$weight_av),
               color = "#e41a1c") +
  geom_segment(aes(x=55+0.2, xend=58+0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 55),]$weight_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 58),]$weight_av),
               color = "#e41a1c") +
  
  annotate("text", x = 7.5, y = 0.5, label = "2021", size = 5) +
  annotate("text", x = 50, y = 0.5, label = "2022", size = 5) +
  
  labs(x = "time (weeks)",
       y = "weight per frond (mg)")

 
px2 <- ggplot(mdata.plot, aes(x=timepoint, y=area_av, color=treatment)) +
  theme_bw(base_size = 14) +
  geom_pointrange(aes(ymin=area_av-area_se, ymax=area_av+area_se), position = position_dodge(0.9)) +
  theme(axis.title.x = element_text(color="black"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))  +
  scale_color_manual(name = "Legend", values=c("#252525", "#4daf4a", "#e41a1c"), labels = c("starting", "Control", "Herbivory"), breaks = c("Starting", "Control", "Herbivory")) +
  scale_x_continuous(breaks = seq(0, 60, by = 4)) +
  expand_limits(x = c(0,61), y = 0) +
  
  geom_rect(data=NULL,aes(xmin = 16, xmax = 40 ,ymin=-Inf,ymax=Inf), fill="gray", color = "transparent") +
  
  geom_segment(aes(x=4-0.2, xend=8-0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 4),]$area_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 8),]$area_av),
               color = "#4daf4a") +
  geom_segment(aes(x=8-0.2, xend=12-0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 8),]$area_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 12),]$area_av),
               color = "#4daf4a") +
  
  geom_segment(aes(x=4+0.2, xend=8+0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 4),]$area_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 8),]$area_av),
               color = "#e41a1c") +
  geom_segment(aes(x=8+0.2, xend=12+0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 8),]$area_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 12),]$area_av),
               color = "#e41a1c") +
  
  geom_segment(aes(x=42-0.2, xend=46-0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 42),]$area_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 46),]$area_av),
               color = "#4daf4a") +
  geom_segment(aes(x=46-0.2, xend=50-0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 46),]$area_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 50),]$area_av),
               color = "#4daf4a") +
  geom_segment(aes(x=50-0.2, xend=55-0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 50),]$area_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 55),]$area_av),
               color = "#4daf4a") +
  geom_segment(aes(x=55-0.2, xend=58-0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 55),]$area_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Control" & mdata.plot$timepoint == 58),]$area_av),
               color = "#4daf4a") +
  
  geom_segment(aes(x=42+0.2, xend=46+0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 42),]$area_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 46),]$area_av),
               color = "#e41a1c") +
  geom_segment(aes(x=46+0.2, xend=50+0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 46),]$area_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 50),]$area_av),
               color = "#e41a1c") +
  geom_segment(aes(x=50+0.2, xend=55+0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 50),]$area_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 55),]$area_av),
               color = "#e41a1c") +
  geom_segment(aes(x=55+0.2, xend=58+0.2,
                   y = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 55),]$area_av,
                   yend = mdata.plot[which(mdata.plot$treatment == "Herbivory" & mdata.plot$timepoint == 58),]$area_av),
               color = "#e41a1c") +
  
  annotate("text", x = 7.5, y = 0.5, label = "2021", size = 5) +
  annotate("text", x = 50, y = 0.5, label = "2022", size = 5) +
  
  labs(x = "time (weeks)",
       y = "area per frond (cm2)")

px <- ggarrange(px2, px1, 
                nrow = 2, labels = c("a.", "b."))

px
```

## Table S2

```R
df <- df[which(df$timepoint != 0),]

model <- lmer(area ~ treatment * timepoint * (1|pond), data = df)
Anova(model)

model <- lmer(weight_fr ~ treatment * timepoint * (1|pond), data = df)
Anova(model)
```

## Main text results

```R
model <- lmer(area ~ treatment * (1|pond) * (1|timepoint), data = df)
Anova(model)

model <- lmer(weight_fr ~ treatment * (1|pond) * (1|timepoint), data = df)
Anova(model)
```

