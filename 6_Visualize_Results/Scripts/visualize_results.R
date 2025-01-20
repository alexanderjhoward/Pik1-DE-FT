# Set up
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
library(tidyverse)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggmsa")
library(ggmsa)
library(ggpubr)
library(caret)
library(stringr)
library(RColorBrewer)

#From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

### Figure 1: Fine-tuned ESM-2 model performance on validation data

# Load validation predictions
pikc_pred <- read.csv('../../2_Model_DE_Data/Output/AvrPikC/esm2_t6_8M_UR50D_regression_predictions.csv')
pikf_pred <- read.csv('../../2_Model_DE_Data/Output/AvrPikF/esm2_t6_8M_UR50D_regression_predictions.csv')
pikc_pred <- pikc_pred %>% 
  mutate(Binding = ifelse(enrichment > 0, "Enriched", "Depleted")) %>%
  mutate(Binding = factor(Binding, levels = c("Enriched", "Depleted")))
pikf_pred <- pikf_pred %>% 
  mutate("Binding" = ifelse(enrichment > 0, "Enriched", "Depleted")) %>%
  mutate(Binding = factor(Binding, levels = c("Enriched", "Depleted")))

# Scatterplot of Avr-PikC predictions
pikc_pred_plot <- ggplot(pikc_pred, aes(enrichment, Predicted_Enrichment)) + 
  geom_smooth(method='lm', formula= y~x, color='black', lwd=0.5) +
  geom_point(aes(fill=Binding), shape=21, color='black', size=1.15, stroke=0.25) +
  stat_cor(method="spearman", label.x.npc = 0.45, label.y.npc = 0.03, size=2.35) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  xlab("Enrichment score (ES)") +
  ylab("Predicted enrichment score (pES)") +
  scale_fill_manual(values = c("#56B4E9", "#FC8D62")) +
  labs(fill="Sequence enrichment\nunder 1 \U03BCM Avr-PikC\nselection") +
  theme_bw() + 
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),  
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"))
ggsave('../Output/Fig1c.png', pikc_pred_plot, height = 56, width = 87, units="mm")

# Scatterplot of Avr-PikF predictions
pikf_pred_plot <- ggplot(pikf_pred, aes(enrichment, Predicted_Enrichment)) + 
  geom_smooth(method='lm', formula= y~x, color='black', lwd=0.5) +
  geom_point(aes(fill=Binding), shape=21, color='black', size=1.15, stroke=0.25) +
  stat_cor(method="spearman", label.x.npc = 0.45, label.y.npc = 0.03, size=2.35) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  xlab("Enrichment score (ES)") +
  ylab("Predicted enrichment score (pES)") +
  scale_fill_manual(values = c("#56B4E9", "#FC8D62")) +
  labs(fill="Sequence enrichment\nunder 1 \U03BCM Avr-PikF\nselection") +
  theme_bw() + 
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),  
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"))
ggsave('../Output/Fig1d.png', pikf_pred_plot, height = 56, width = 87, units="mm")



### Figure 2: MSA of sequences detected in the 3k RGP dataset

# Multiple sequence alignment
msa <- ggmsa(paste0('../../3_Variant_Calling/Output/Pikp1_HMA_Variants_Unique_Aligned_Formatted.fa'), 1, 82, seq_name=TRUE, char_width = 0.55, color = "Zappo_AA", border = "white", consensus_views = TRUE, use_dot = TRUE, ref = "Pikh-1_(24)") + 
  scale_x_continuous(breaks = seq(0, 82, by = 10), expand = c(0, 0)) +
  scale_y_discrete(labels=function(x) gsub("_", " ", x, fixed=TRUE)) +
  theme(text=element_text(size=5.25))
ggsave(paste0('../Output/Fig2a.png'), msa, height = 31.5, width = 158, units="mm", dpi=1000)

# Load predicted variant enrichment scores
seq_order <- c("Pikh-1", "Pikp-1", "Bangla", "Larha_Mugad", "Chundi", 
               "Periya_Vellai", "Vellai_Kolomban", "Taipei_Woo_Co", "Adukkan", 
               "Xintuan_Hei_Gu", "Zimangfeie", "Banjiemang2", "Sanhuangzhan_No_2")
pES_C <- read.csv('../../4_G2P_Predictions/Output/AvrPikC_esm2_t6_8M_UR50D_predictions.csv') %>%
  mutate(Ligand = 'Avr-PikC') %>%
  mutate(pES = Predicted_AvrPikC_Enrichment) %>%
  select(Variety, Ligand, pES)
pES_F <- read.csv('../../4_G2P_Predictions/Output/AvrPikF_esm2_t6_8M_UR50D_predictions.csv') %>%
  mutate(Ligand = 'Avr-PikF') %>%
  mutate(pES = Predicted_AvrPikF_Enrichment) %>%
  select(Variety, Ligand, pES)
pES_heatmap <- rbind(pES_C, pES_F) %>%
  mutate(Variety = factor(Variety, levels = rev(seq_order)))
rm(pES_C, pES_F, seq_order)

# Plot pES values as a heatmap
pES_heatmap_plot <- ggplot(pES_heatmap, aes(Ligand, Variety, fill= pES)) + 
  geom_tile(color='white') +
  geom_text(aes(label=sprintf("%0.2f", pES)), color = "black", size = 1.45) +
  scale_fill_gradient2(low = "#FC8D62", midpoint = 0, mid = "white", high = "#56B4E9") +
  scale_x_discrete(expand = c(0.01, 0.01)) +
  scale_y_discrete(expand = c(0.01, 0.01)) +
  coord_cartesian(expand = F) +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=4.5),
        plot.margin=unit(c(1,1,9,4),"mm"),
        legend.position = c(0.39,-0.3), 
        legend.direction = "horizontal",
        legend.margin=margin(t = 0, unit='cm'),
        legend.key.width = unit(0.65,"line"),
        legend.key.height = unit(0.35,"line"),
        legend.text = element_text(size=4.5), 
        legend.title = element_text(size=4.5))
ggsave('../Output/Fig2b.png', pES_heatmap_plot, height = 35, width = 25, units="mm", dpi=1000)



### Figure 4: Fine-tuned ESM-2 model performance on NUDT15 data

# Plot predictions on validation data
preds <- read.csv('../../5_NUDT15_Modeling/Output/esm2_t6_8M_UR50D_regression_predictions.csv') %>% 
  mutate(Final_classification = ifelse(Score > 0, "Functional", "Nonfunctional"))
pred_plot <- ggplot(preds, aes(Score, Predicted_Score)) + 
  geom_smooth(method='lm', formula= y~x, color='black', lwd=0.5) +
  geom_point(aes(fill=Final_classification), shape=21, color='black', size=0.95, stroke=0.25) +
  stat_cor(method="spearman", label.x.npc = 0.45, label.y.npc = 0.02, size=2.35) +
  scale_x_continuous(breaks = c(-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6), expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4), expand = c(0.01, 0.01)) +
  xlab("Functionality score (FS)") +
  ylab("Predicted functionality score (pFS)") +
  scale_fill_manual(values = c("#56B4E9","#FC8D62")) + 
  labs(fill="NUDT15 functionality") +
  theme_bw() + 
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),  
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"))
ggsave('../Output/Fig4a.png', pred_plot, height = 58, width = 90, units="mm")

# Plot predictions for known clinical variants
## Merge predictions made by modeling with scores from paper
var_preds <- read.csv('../../5_NUDT15_Modeling/Output/esm2_t6_8M_UR50D_regression_variant_predictions.csv') %>%
  select(variant, seq, effect, Predicted_Score, Allele.Count, Allele.Frequency)
scores <- read.csv('../../5_NUDT15_Modeling/Source/NUDT15.csv') %>%
  drop_na(Final.NUDT15.activity.Score) %>%
  mutate(Score=Final.NUDT15.activity.Score) %>%
  select(Protein, Score, Final_classification) %>% # After calculating scores, only keep sequences that have more than n reads in the starting or ending library
  mutate(Score = log(Score + 0.55)) # Minimum cut-off for wt-like activity set by authors was 0.45, so add 0.55 to make all wt-like activity >= 1
var_preds <- merge(var_preds, scores, by.x="seq", by.y="Protein", all=T) %>%
  drop_na(effect)
## Format data for heat map
clin_vars_true <- var_preds %>%
  filter(effect %in% c("Benign", "Toxic")) %>%
  mutate(variant = gsub("WT", "Wild type", variant)) %>%
  arrange(effect) %>%
  select(variant, effect, Score) %>%
  mutate(group = "FS")
clin_vars_pred <- var_preds %>%
  filter(effect %in% c("Benign", "Toxic")) %>%
  mutate(variant = gsub("WT", "Wild type", variant)) %>%
  arrange(effect) %>%
  mutate(Score = Predicted_Score) %>%
  select(variant, effect, Score) %>%
  mutate(group = "pFS")
clin_vars_heatmap <- rbind(clin_vars_true, clin_vars_pred) %>%
  mutate(variant = factor(variant, levels = rev(c('Wild type', 'Q6E', 'V18I', 'V93I', 'R11Q', 'S83Y', 
                                                  'K33E', 'R139H', 'R139C', 'R34T', 'V75G', 'G17V18del', 
                                                  'G17V18dup', 'G17V18dup/R139C')))) %>%
  mutate(group = factor(group, levels = c('FS', 'pFS')))
rm(clin_vars_pred, clin_vars_true)
## Plot a heatmap comparing function scoring and model scoring of each clinical variant
clin_vars_heatmap_plot <- ggplot(clin_vars_heatmap, aes(group, variant, fill= Score)) + 
  geom_tile(color='white') +
  geom_text(aes(label=sprintf("%0.2f", Score)), color = "black", size = 2.5) +
  scale_fill_gradient2(low = "#FC8D62", midpoint = 0, mid = "white", high = "#56B4E9") +
  scale_x_discrete(expand = c(0.01, 0.01)) +
  scale_y_discrete(expand = c(0.01, 0.01)) +
  coord_cartesian(expand = F) +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_text(size=7),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=7),
        plot.margin=unit(c(1,1,9,4),"mm"),
        legend.position = c(0.33,-0.2), 
        legend.direction = "horizontal",
        legend.margin=margin(t = 0, unit='cm'),
        legend.key.width = unit(1.15,"line"),
        legend.key.height = unit(0.35,"line"),
        legend.text = element_text(size=7), 
        legend.title = element_text(size=7))
ggsave('../Output/Fig4b.png', clin_vars_heatmap_plot, height = 58, width = 60, units="mm", dpi=1000)

# Plot model scoring of variants with no known clinical characterization or function score
unchar_preds <- var_preds %>%
  filter(effect=="Uncharacterized") %>%
  filter(!(Final_classification %in% c("wt-like", "damaging"))) %>%
  separate(variant, into = c("wt_res", "intermediate"), sep = "(?<=[A-Za-z])(?=[0-9])", remove=F) %>%
  separate(intermediate, into = c("pos", "extra"), sep = "(?<=[0-9])(?=[A-Za-z])") %>%
  mutate(pos = as.integer(pos)) %>%
  mutate(mut_type = ifelse(str_detect(variant, "del"), "Deletion",
                           ifelse(str_detect(variant, "ins"), "Insertion",
                                  ifelse(str_detect(variant, "dup"), "Duplication", "Substitution")))) %>%
  select(seq, variant, pos, Predicted_Score, mut_type, Allele.Count, Allele.Frequency)

pheno_plot <- ggplot(unchar_preds, aes(pos, Predicted_Score, ymax = Predicted_Score, ymin = 0)) + 
  geom_pointrange(aes(fill=mut_type, size=Allele.Count), shape=21, color='black', stroke=0.25, lwd=0.25, linetype = 2) +
  scale_fill_brewer(palette = 'Set2') +
  geom_segment(aes(x = 1, y = 0, xend = 164, yend = 0), lwd=0.75) +
  scale_x_continuous(breaks = c(1, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 164), expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.05, 0.05)) +
  scale_size(range = c(0.35,1.35), breaks = c(1, 5, 20)) +
  labs(x = "NUDT15 Position", y = "Predicted functionality score (pFS)", fill="Mutation type", size="Allele count") +
  theme_bw() +
  theme(axis.title.y=element_text(size=8),
        axis.text.y=element_text(size=6),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=6),
        legend.text = element_text(size=6), 
        legend.title = element_text(size=8),
        legend.key.height= unit(3.75, 'mm'),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(2, "mm"),
        legend.margin = margin(0, 0, 0, 0),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25))
ggsave('../Output/Fig4c.png', pheno_plot, height = 51, width = 150, units="mm", dpi=1000)
