# Set up
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
library(tidyverse)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggmsa")
library(ggmsa)

### Figure 1: Fine-tuned ESM-2 model performance on validation data

# Load validation predictions
pikc_pred <- read.csv('../../2_Model_DE_Data/Output/AvrPikC/esm2_t6_8M_UR50D_regression_predictions.csv')
pikf_pred <- read.csv('../../2_Model_DE_Data/Output/AvrPikF/esm2_t6_8M_UR50D_regression_predictions.csv')
pikc_pred <- pikc_pred %>% 
  mutate("Binding" = ifelse(enrichment > 0, "Positive", "Negative"))
pikf_pred <- pikf_pred %>% 
  mutate("Binding" = ifelse(enrichment > 0, "Positive", "Negative"))

# Scatterplot of Avr-PikC predictions
pikc_pred_plot <- ggplot(pikc_pred, aes(enrichment, Predicted_Enrichment)) + 
  geom_smooth(method='lm', formula= y~x, color='black', lwd=0.5) +
  geom_point(aes(fill=Binding), shape=21, color='black', size=1.15, stroke=0.25) +
  stat_cor(method="spearman", label.x.npc = 0.45, label.y.npc = 0.03, size=2.35) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  xlab("Enrichment score (ES)") +
  ylab("Predicted enrichment score (pES)") +
  scale_fill_manual(values = c("#FC8D62","#66C2A5")) + 
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
  scale_fill_manual(values = c("#FC8D62","#66C2A5")) + 
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
  scale_fill_gradient2(low = "#FC8D62", midpoint = 0, mid = "white", high = "#66C2A5") +
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
