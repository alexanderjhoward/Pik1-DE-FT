# Set up
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
library(tidyverse)
library(seqinr)
library(ggpubr)
library(caret)
library(stringr)
library(RColorBrewer)

# Load and save  wt information
wt <- read.csv(file='NUDT15.csv') %>%
  filter(Info=="WT")
wtseq <- wt$Protein[1]
rm(wt)

# List of variants with known clinical phenotypes
## Toxic variants according to literature search
# Source 1: https://www.nature.com/articles/ng.3508
# Source 2: https://jamanetwork.com/journals/jama/fullarticle/2725687
# Source 3: https://pmc.ncbi.nlm.nih.gov/articles/PMC7071893/#sec11
toxic <- c('Arg139Cys', 'Arg139His', 'Lys33Glu', 'Arg34Thr', 'Val75Gly')
toxic_df <- as.data.frame(toxic) %>%
  separate(toxic, into = c("wt_res", "intermediate"), sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  separate(intermediate, into = c("pos", "mut_res"), sep = "(?<=[0-9])(?=[A-Za-z])") %>%
  mutate(pos = as.integer(pos)) %>%
  mutate(wt_res = a(wt_res)) %>%
  mutate(mut_res = a(mut_res)) %>%
  mutate(seq = wtseq) %>%
  mutate(effect = "Toxic")
for(i in 1:nrow(toxic_df)){
  substring(toxic_df[i,4], (toxic_df[i,2]), (toxic_df[i,2])) <- toxic_df[i,3]
}
## Benign variants according to literature search
benign <- c('Val18Ile', 'Gln6Glu', 'Arg11Gln', 'Ser83Tyr', 'Val93Ile')
benign_df <- as.data.frame(benign) %>%
  separate(benign, into = c("wt_res", "intermediate"), sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  separate(intermediate, into = c("pos", "mut_res"), sep = "(?<=[0-9])(?=[A-Za-z])") %>%
  mutate(pos = as.integer(pos)) %>%
  mutate(wt_res = a(wt_res)) %>%
  mutate(mut_res = a(mut_res)) %>%
  mutate(seq = wtseq) %>%
  mutate(effect = "Benign")
for(i in 1:nrow(benign_df)){
  substring(benign_df[i,4], (benign_df[i,2]), (benign_df[i,2])) <- benign_df[i,3]
}
#Gly17Val18del
toxic_df <- rbind(toxic_df, data.frame(wt_res="GV", pos=17, mut_res="", seq=paste0(substring(wtseq, 1, 16), substring(wtseq, 19)), effect="Toxic"))
#Gly17Val18dup
toxic_df <- rbind(toxic_df, data.frame(wt_res="GV", pos=17, mut_res="GVGV", seq=paste0(substring(wtseq, 1, 18), substring(wtseq, 17, 18), substring(wtseq, 19)), effect="Toxic"))
#Gly17Val18_dup/Arg139Cys
toxic_df <- rbind(toxic_df, data.frame(wt_res="GV/A", pos=17, mut_res="GV/C", seq=paste0(substring(wtseq, 1, 18), substring(wtseq, 17, 18), substring(wtseq, 19)), effect="Toxic"))
substring(toxic_df[8,4], 139+2, 139+2) <- "C"
#WT
benign_df <- rbind(benign_df, data.frame(wt_res="", pos=0, mut_res="", seq=wtseq, effect="Benign"))
# Collect allele count info from genomAD
genomad <- read.csv('gnomAD_v4.1.0_ENSG00000136159_2025_01_13_14_21_37.csv') %>%
  mutate(Protein.Consequence = gsub('^p.', '', Protein.Consequence))
genomadvars <- subset(genomad, genomad$Protein.Consequence %in% c(toxic, benign)) %>%
  select(Protein.Consequence, Allele.Count, Allele.Frequency) %>%
  separate(Protein.Consequence, into = c("wt_res", "intermediate"), sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  separate(intermediate, into = c("pos", "mut_res"), sep = "(?<=[0-9])(?=[A-Za-z])") %>%
  mutate(pos = as.integer(pos)) %>%
  mutate(wt_res = a(wt_res)) %>%
  mutate(mut_res = a(mut_res)) %>%
  mutate(variant = paste0(wt_res, pos, mut_res)) %>%
  select(variant, Allele.Count, Allele.Frequency)
genomadvars_extras <- subset(genomad, genomad$Protein.Consequence %in% c('Gly17_Val18del', 'Gly17_Val18dup')) %>%
  select(Protein.Consequence, Allele.Count, Allele.Frequency) %>%
  mutate(variant = Protein.Consequence) %>%
  mutate(variant = gsub("Gly17_Val18del", "G17V18del", variant)) %>%
  mutate(variant = gsub("Gly17_Val18dup", "G17V18dup", variant)) %>%
  select(variant, Allele.Count, Allele.Frequency)
genomadvars_missing <- data.frame(variant = c("WT", "G17V18dup/R139C"), 
                                  Allele.Count = c(NA, NA), 
                                  Allele.Frequency = c(NA, NA))
genomadvars_total <- rbind(genomadvars, genomadvars_extras, genomadvars_missing)
rm(genomadvars, genomadvars_extras, genomadvars_missing)
## Collect all variants
variants <- rbind(benign_df, toxic_df) %>%
  mutate(variant = paste0(wt_res, pos, mut_res)) %>%
  mutate(variant = gsub("0","WT",variant)) %>%
  mutate(variant = gsub("GV17GVGV","G17V18dup",variant)) %>%
  mutate(variant = gsub("GV17","G17V18del",variant)) %>%
  mutate(variant = gsub("GV/A17GV/C","G17V18dup/R139C",variant)) %>%
  select(variant, seq, effect) %>%
  merge(genomadvars_total, by="variant") %>%
  arrange(effect)
rm(benign_df, toxic_df, benign, toxic, genomadvars_total)


# List of variants without a known phenotype
## Missense variants
genomad_missense <- genomad %>%
  select(Protein.Consequence, VEP.Annotation, ClinVar.Germline.Classification, Allele.Count, Allele.Frequency) %>%
  filter(VEP.Annotation == 'missense_variant') %>%
  mutate(Protein.Consequence = gsub('^p.', '', Protein.Consequence)) %>%
  separate(Protein.Consequence, into = c("wt_res", "intermediate"), sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  separate(intermediate, into = c("pos", "mut_res"), sep = "(?<=[0-9])(?=[A-Za-z])") %>%
  mutate(wt_res = gsub("p.", "", wt_res)) %>%
  mutate(wt_res = a(wt_res)) %>%
  mutate(mut_res = a(mut_res)) %>%
  mutate(pos = as.integer(pos)) %>%
  mutate(seq = wtseq) %>%
  mutate(effect = "Uncharacterized") %>%
  select(wt_res, pos, mut_res, seq, effect, Allele.Count, Allele.Frequency) 
for(i in 1:nrow(genomad_missense)){
  substring(genomad_missense[i,4], (genomad_missense[i,2]), (genomad_missense[i,2])) <- genomad_missense[i,3]
}
genomad_missense <- genomad_missense %>%
  mutate(variant = paste0(wt_res, pos, mut_res)) %>%
  filter(!(seq %in% variants$seq)) %>%
  select(variant, seq, effect, Allele.Count, Allele.Frequency) %>%
  distinct(seq, .keep_all = TRUE)
## Single deletions
genomad_dels <- genomad %>%
  select(Protein.Consequence, VEP.Annotation, ClinVar.Germline.Classification, Allele.Count, Allele.Frequency) %>%
  filter(VEP.Annotation %in% c('inframe_deletion', 'inframe_insertion')) %>%
  filter(!(grepl("_", Protein.Consequence))) %>%
  mutate(Protein.Consequence = gsub('^p.', '', Protein.Consequence)) %>%
  separate(Protein.Consequence, into = c("wt_res", "intermediate"), sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  separate(intermediate, into = c("pos", "mut_res"), sep = "(?<=[0-9])(?=[A-Za-z])") %>%
  mutate(wt_res = gsub("p.", "", wt_res)) %>%
  mutate(wt_res = a(wt_res)) %>%
  mutate(pos = as.integer(pos)) %>%
  mutate(seq = wtseq) %>%
  mutate(effect = "Uncharacterized") %>%
  select(wt_res, pos, mut_res, seq, effect, Allele.Count, Allele.Frequency) 
for(i in 1:nrow(genomad_dels)){
  genomad_dels[i,4] <- paste0(substring(wtseq, 1, genomad_dels[i,2]-1), substring(wtseq, genomad_dels[i,2]+1))
}
genomad_dels <- genomad_dels %>%
  mutate(variant = paste0(wt_res, pos, mut_res)) %>%
  filter(!(seq %in% variants$seq)) %>%
  select(variant, seq, effect, Allele.Count, Allele.Frequency) %>%
  distinct(seq, .keep_all = TRUE)
## Range mutations
genomad_range <- genomad %>%
  select(Protein.Consequence, VEP.Annotation, ClinVar.Germline.Classification, Allele.Count, Allele.Frequency) %>%
  filter(VEP.Annotation %in% c('inframe_deletion', 'inframe_insertion')) %>%
  filter(grepl("_", Protein.Consequence)) %>%
  mutate(Protein.Consequence = gsub('^p.', '', Protein.Consequence)) %>%
  separate(Protein.Consequence, c('First', 'Second'), sep="_") %>%
  separate(First, into = c("first_wt_res", "first_pos"), sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  separate(Second, into = c("second_wt_res", "intermediate"), sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  separate(intermediate, into = c("second_pos", "mut_res"), sep = "(?<=[0-9])(?=[A-Za-z])") %>%
  mutate(first_wt_res = a(first_wt_res)) %>%
  mutate(second_wt_res = a(second_wt_res)) %>%
  mutate(first_pos = as.integer(first_pos)) %>%
  mutate(second_pos = as.integer(second_pos)) %>%
  mutate(seq = wtseq) %>%
  mutate(effect = "Uncharacterized") %>%
  select(first_wt_res, first_pos, second_wt_res, second_pos, mut_res, seq, effect, Allele.Count, Allele.Frequency)
genomad_range_dup <- genomad_range %>%
  filter(mut_res=='dup')
for(i in 1:nrow(genomad_range_dup)){
  genomad_range_dup[i,6] <- paste0(substring(wtseq, 1, genomad_range_dup[i,4]), substring(wtseq, genomad_range_dup[i,2], genomad_range_dup[i,4]), substring(wtseq, genomad_range_dup[i,4]+1))
}
genomad_range_dup <- genomad_range_dup %>%
  mutate(variant = paste0(first_wt_res, first_pos, second_wt_res, second_pos, mut_res)) %>%
  filter(!(seq %in% variants$seq)) %>%
  select(variant, seq, effect, Allele.Count, Allele.Frequency) %>%
  distinct(seq, .keep_all = TRUE)
genomad_range_del <- genomad_range %>%
  filter(mut_res=='del')
for(i in 1:nrow(genomad_range_del)){
  genomad_range_del[i,6] <- paste0(substring(wtseq, 1, genomad_range_del[i,2]-1), substring(wtseq, genomad_range_del[i,4]+1))
}
genomad_range_del <- genomad_range_del %>%
  mutate(variant = paste0(first_wt_res, first_pos, second_wt_res, second_pos, mut_res)) %>%
  filter(!(seq %in% variants$seq)) %>%
  select(variant, seq, effect, Allele.Count, Allele.Frequency) %>%
  distinct(seq, .keep_all = TRUE)
genomad_ins <- genomad_range %>%
  filter(!(mut_res %in% c('del', 'dup')))
genomad_ins[1,6] <- paste0(substring(wtseq, 1, genomad_ins[1,2]), "RV", substring(wtseq, genomad_ins[1,4]))
genomad_ins <- genomad_ins %>%
  mutate(variant = paste0(first_wt_res, first_pos, second_wt_res, second_pos, "insRV")) %>%
  filter(!(seq %in% variants$seq)) %>%
  select(variant, seq, effect, Allele.Count, Allele.Frequency)
## Collect variants
genomad_vars <- rbind(genomad_missense, genomad_dels, genomad_ins, genomad_range_del, genomad_range_dup) %>%
  arrange(effect)
rm(genomad_missense, genomad_dels, genomad_ins, genomad_range_del, genomad_range_dup, genomad_range, genomad)

# Load and clean data
# Adjust scoring such that all wild-type sequences are positive values and all LOF sequences are negative values
df <- read.csv(file='NUDT15.csv') %>%
  drop_na(Final.NUDT15.activity.Score) %>%
  mutate(Score=Final.NUDT15.activity.Score) %>%
  select(Protein, Score, Final_classification) %>% # After calculating scores, only keep sequences that have more than n reads in the starting or ending library
  mutate(Score = log(Score + 0.55)) %>% # Minimum cut-off for wt-like activity set by authors was 0.45, so add 0.55 to make all wt-like activity >= 1
  filter(!(Protein %in% variants$seq)) %>%
  mutate(bin = ntile(Score, round(n()/100, digits = 0)))
write.csv(df, file='NUDT15_cleaned_adj.csv', row.names = F)

# Remove any genomAD variants that have been already tested in this assay
genomad_vars <- genomad_vars %>%
  filter(!(seq %in% df$Protein))
full_vars <- rbind(variants, genomad_vars)
rm(variants, genomad_vars)
write.csv(full_vars, file='NUDT15_variants.csv', row.names = F)


#################
# After modeling#
#################


# Plot predictions on validation data
preds <- read.csv('esm2_t6_8M_UR50D_regression_predictions.csv') %>% 
  mutate(Final_classification = ifelse(Score > 0, "Functional", "Nonfunctional"))
pred_plot <- ggplot(preds, aes(Score, Predicted_Score)) + 
  geom_smooth(method='lm', formula= y~x, color='black', lwd=0.5) +
  geom_point(aes(fill=Final_classification), shape=21, color='black', size=0.95, stroke=0.25) +
  stat_cor(method="spearman", label.x.npc = 0.45, label.y.npc = 0.02, size=2.35) +
  scale_x_continuous(breaks = c(-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6), expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4), expand = c(0.01, 0.01)) +
  xlab("Functionality score (FS)") +
  ylab("Predicted functionality score (pFS)") +
  scale_fill_manual(values = c("#66C2A5","#FC8D62")) + 
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
ggsave('Fig4a.png', pred_plot, height = 58, width = 90, units="mm")

# Plot predictions made for known clinical variants
## Merge predictions made by modeling with scores from paper
var_preds <- read.csv('esm2_t6_8M_UR50D_regression_variant_predictions.csv') %>%
  select(variant, seq, effect, Predicted_Score, Allele.Count, Allele.Frequency)
scores <- read.csv('NUDT15.csv') %>%
  drop_na(Final.NUDT15.activity.Score) %>%
  mutate(Score=Final.NUDT15.activity.Score) %>%
  select(Protein, Score, Final_classification) %>% # After calculating scores, only keep sequences that have more than n reads in the starting or ending library
  mutate(Score = log(Score + 0.55)) # Minimum cut-off for wt-like activity set by authors was 0.45, so add 0.55 to make all wt-like activity >= 1
var_preds <- merge(var_preds, scores, by.x="seq", by.y="Protein", all=T) %>%
  drop_na(effect)
## Identify the clinical variants to plot
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
  scale_fill_gradient2(low = "#FC8D62", midpoint = 0, mid = "white", high = "#66C2A5") +
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
ggsave('Fig4b.png', clin_vars_heatmap_plot, height = 58, width = 60, units="mm", dpi=1000)

# Plot model scoring of variants with no known characterization or function score
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
ggsave('Fig4c.png', pheno_plot, height = 51, width = 150, units="mm", dpi=1000)
