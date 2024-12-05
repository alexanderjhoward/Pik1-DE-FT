# Set up
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(seqinr)
library(ggdendro)
library(dendextend)
library(ggmsa)

### Gather all unique HMA sequence variants
# Format and combine searched Pik-1 HMA domain exons
format_pik_hma <- function(a, b){
  # Load in sample annotations
  ids <- read.csv('../Source/accession_info.csv')
  annots <- read.table('../Source/3K_list_sra_ids.txt', sep='\t', header=T) %>%
    mutate(Country_Origin_updated = gsub(" ", "'", Country_Origin_updated)) %>%
    mutate(Country_Origin_updated = gsub("Korea,_Republic_of", "Republic_of_Korea", Country_Origin_updated)) %>%
    mutate(Country_Origin_updated = gsub("_no_info", "Not_Available", Country_Origin_updated))
  sample_info <- merge(annots, ids, by.x="SRA.Accession", by.y="SRA")
  # Load in and format exons
  exon1 <- read.csv(paste0('../../3_Search_3KRGP/Output/', a, '/', a, '_E1_pep.csv'))
  exon1 <- separate(exon1, id, into=c('BIOSAMPLE', "Exon1"), sep="_")
  exon2 <- read.csv(paste0('../../3_Search_3KRGP/Output/', a, '/', a, '_E2_pep.csv'))
  exon2 <- separate(exon2, id, into=c('BIOSAMPLE', "Exon2"), sep="_")
  # Merge exons into contiguous sequence
  pik <- merge(sample_info, exon1, by='BIOSAMPLE')
  pik <- merge(pik, exon2, by='BIOSAMPLE')
  # Filter out HMA domains called from heterozygous samples and format data
  het_reads <- c('GLKQKIVFKIPMEGNNCRSKAMALVASTGGVDSVALVGDLRDKIEVVGDGIDSINLVSALRKKVGPAMFLEVSQAKED',
                 'GLKQTIVIKVAMEGNNCRSKAMDLVKSTGGVYSVSLAGDLRDKIEVVGYGIDPIKLISALRKKVGHAELLQVSQAKKD',
                 'GGEMQKIVFKIPMVDDKSRTKAMSLVASTVGVHSVAIAGDLRDEVVVVGYGIDPINLVSALRKKVGPAMFLEVSQAKED',
                 'GGEMQKIVFKIPMVDDKSRTKAMSLVASTVGVHSVAIAGDLRDDVVVVGDGIDSINLVSALRKKVDPAMFLEVSQAKED')
  pik <- pik %>%
    mutate(SEQ = paste(seq.x, seq.y, sep="")) %>%
    filter(!(SEQ %in% het_reads)) %>%
    select(c("BIOSAMPLE", "SRA.Accession", "X3K_DNA_IRIS_UNIQUE_ID", "Genetic_Stock_varname", "Country_Origin_updated", "SEQ"))
  colnames(pik) <- c("BIOSAMPLE", "SRA", "IRIS", "STOCK", "COUNTRY", "SEQ")
  # Back-fill the C-terminal deletion found in 2 samples
  pik <- pik %>% mutate(SEQ = ifelse(SEQ == "GLKQKIVIKVAMEGNNCRSKAMALVASTGGVDSVALVGDLRDKIVVVGYGIDPIKLISALRKKVGDAELLQVSD", paste0(SEQ, "VKET"), SEQ))
  # Save sequence results for all searched varieties
  write.csv(pik, file = paste0('../Output/Data/', b, '_HMA_Variants.csv'),  row.names = FALSE)
  # Identify frequency of variants
  pik_freq <- as.data.frame(table(pik$SEQ))
  # Filter for unique variants
  pik_uniq <- pik %>%
    filter(!grepl("X", SEQ)) %>%
    filter(!grepl("\\*", SEQ)) %>%
    distinct(SEQ, .keep_all = TRUE)
  pik_uniq <- merge(pik_uniq, pik_freq, by.x="SEQ", by.y="Var1") %>% 
    separate(STOCK, c("Variety", "ID"), sep="::") %>%
    mutate(Variety = gsub("(?<=\\p{L})(\\p{L}+)", "\\L\\1", Variety, perl = TRUE))
  return(pik_uniq)
}
pikm1_uniq <- format_pik_hma('Kit', 'Pikm1')
pikp1_uniq <- format_pik_hma('N22', 'Pikp1')

# Format names of variants and save
## Pikm-1 variants
pikm1_uniq <- pikm1_uniq  %>%
  mutate(Variety = gsub("Pdr", "PDR", Variety)) %>%
  mutate(Variety = gsub("Wir", "WIR", Variety)) %>%
  mutate(Variety = gsub("Uprh", "UPRH", Variety)) %>%
  mutate(Variety = gsub("Inia", "INIA", Variety)) %>%
  mutate(Variety = gsub("81_A_32", "Pik*-1", Variety)) %>%
  mutate(Variety = gsub("\\(Nam\\)", "Nam", Variety)) %>%
  mutate(Variety = gsub("INIA_Tacuari", "Pikm-1", Variety)) %>%
  mutate(Variety = gsub("Sican", "Piks-1", Variety)) 
write.csv(pikm1_uniq, file = '../Output/Data/Pikm1_HMA_Variants_Unique.csv',  row.names = FALSE)
## Pikp-1 variants
pikp1_uniq <- pikp1_uniq %>%
  mutate(Variety = gsub("Iet_14720", "Pikh-1", Variety)) %>%
  mutate(Variety = gsub("Simul_Khuri", "Pikp-1", Variety))
write.csv(pikp1_uniq, file = '../Output/Data/Pikp1_HMA_Variants_Unique.csv',  row.names = FALSE)

# Save FASTA of unique variants in specified order (determined by a dendogram that's generated further down)
write_fasta <- function(a, b, c){
  fasta <- a %>%
    mutate(Variety = factor(Variety, levels = b)) %>%
    arrange(Variety) %>%
    mutate(Number = sprintf("%02d", 1:nrow(a))) %>%
    mutate(Variety = paste0('>', Number, ";", Variety, "_(", Freq, ")")) %>%
    select(Variety, SEQ)
  fasta <- do.call(rbind, lapply(seq(nrow(fasta)), function(i) t(fasta[i, ])))
  write.table(fasta, file = paste0('../Output/Data/', c, '_HMA_Variants_Unique.fa'),  row.names = FALSE, col.names = FALSE, quote = FALSE)
}
pikm_seq_order <- c("Pikm-1", "Jana_Nam", "571", "Avo", "PDR_34-2-1-2", 
                    "WIR_1072", "UPRH_197", "Piks-1", "Samrong_2", "Pik*-1")
pikp_seq_order <- c("Pikh-1", "Pikp-1", "Bangla", "Larha_Mugad", "Chundi", 
                    "Periya_Vellai", "Vellai_Kolomban", "Taipei_Woo_Co", "Adukkan", 
                    "Xintuan_Hei_Gu", "Zimangfeie", "Banjiemang2", "Sanhuangzhan_No_2")
write_fasta(pikm1_uniq, pikm_seq_order, 'Pikm1')
write_fasta(pikp1_uniq, pikp_seq_order, 'Pikp1')

#################################################################
# I aligned and formatted these FASTA files on the command line #
#################################################################

### Gather all unique Pik1/2 variants possessing San/Vel HMA domains
format_pik_vars <- function(){
  # Load in sample annotations
  ids <- read.csv('../Source/accession_info.csv')
  annots <- read.table('../Source/3K_list_sra_ids.txt', sep='\t', header=T) %>%
    mutate(Country_Origin_updated = gsub(" ", "'", Country_Origin_updated)) %>%
    mutate(Country_Origin_updated = gsub("Korea,_Republic_of", "Republic_of_Korea", Country_Origin_updated)) %>%
    mutate(Country_Origin_updated = gsub("_no_info", "Not_Available", Country_Origin_updated))
  sample_info <- merge(annots, ids, by.x="SRA.Accession", by.y="SRA")
  # Load in and format exons
  p1_exon1 <- read.csv('../../3_Search_3KRGP/Output/Vars/Vars_P1E1_pep.csv')
  p1_exon1 <- separate(p1_exon1, id, into=c('BIOSAMPLE', "Exon1"), sep="_")
  p1_exon2 <- read.csv('../../3_Search_3KRGP/Output/Vars/Vars_P1E2_pep.csv')
  p1_exon2 <- separate(p1_exon2, id, into=c('BIOSAMPLE', "Exon2"), sep="_")
  p1_exon3 <- read.csv('../../3_Search_3KRGP/Output/Vars/Vars_P1E3_pep.csv')
  p1_exon3 <- separate(p1_exon3, id, into=c('BIOSAMPLE', "Exon3"), sep="_")
  p2_exon1 <- read.csv('../../3_Search_3KRGP/Output/Vars/Vars_P2E1_pep.csv')
  p2_exon1 <- separate(p2_exon1, id, into=c('BIOSAMPLE', "Exon1"), sep="_")
  p2_exon2 <- read.csv('../../3_Search_3KRGP/Output/Vars/Vars_P2E2_pep.csv')
  p2_exon2 <- separate(p2_exon2, id, into=c('BIOSAMPLE', "Exon2"), sep="_")
  # Merge exons into contiguous sequences
  p1 <- merge(sample_info, p1_exon1, by='BIOSAMPLE')
  p1 <- merge(p1, p1_exon2, by='BIOSAMPLE')
  p1 <- merge(p1, p1_exon3, by='BIOSAMPLE')
  p1$SEQ <- paste(p1$seq.x, p1$seq.y, p1$seq, sep="")
  p2 <- merge(sample_info, p2_exon1, by='BIOSAMPLE')
  p2 <- merge(p2, p2_exon2, by='BIOSAMPLE')
  p2$SEQ <- paste(p2$seq.x, p2$seq.y, sep="")
  # Format data
  p1 <- p1 %>%
    select(c("BIOSAMPLE", "SRA.Accession", "X3K_DNA_IRIS_UNIQUE_ID", "Genetic_Stock_varname", "Country_Origin_updated", "SEQ"))
  colnames(p1) <- c("BIOSAMPLE", "SRA", "IRIS", "STOCK", "COUNTRY", "SEQ")
  p2 <- p2 %>%
    select(c("BIOSAMPLE", "SRA.Accession", "X3K_DNA_IRIS_UNIQUE_ID", "Genetic_Stock_varname", "Country_Origin_updated", "SEQ"))
  colnames(p2) <- c("BIOSAMPLE", "SRA", "IRIS", "STOCK", "COUNTRY", "SEQ")
  # Save sequence results for all searched varieties
  write.csv(p1, file = '../Output/Data/Vars_Pik1_Sequences.csv',  row.names = FALSE)
  write.csv(p2, file = '../Output/Data/Vars_Pik2_Sequences.csv',  row.names = FALSE)
  # Identify frequency of variants
  p1_freq <- as.data.frame(table(p1$SEQ))
  p2_freq <- as.data.frame(table(p2$SEQ))
  # Filter for unique variants
  p1_uniq <- p1 %>% distinct(SEQ, .keep_all = TRUE)
  p2_uniq <- p2 %>% distinct(SEQ, .keep_all = TRUE)
  p1_uniq <- merge(p1_uniq, p1_freq, by.x="SEQ", by.y="Var1") %>% 
    separate(STOCK, c("Variety", "ID"), sep="::") %>%
    mutate(Variety = gsub("(?<=\\p{L})(\\p{L}+)", "\\L\\1", Variety, perl = TRUE))
  p2_uniq <- merge(p2_uniq, p2_freq, by.x="SEQ", by.y="Var1") %>% 
    separate(STOCK, c("Variety", "ID"), sep="::") %>%
    mutate(Variety = gsub("(?<=\\p{L})(\\p{L}+)", "\\L\\1", Variety, perl = TRUE))
  write.csv(p1_uniq, file = '../Output/Data/Vars_Pik1_Sequences_Unique.csv',  row.names = FALSE)
  write.csv(p2_uniq, file = '../Output/Data/Vars_Pik2_Sequences_Unique.csv',  row.names = FALSE)
}
format_pik_vars()

#################################################################
### Figure 1: Fine-tuned ESM-2 model performance on validation data
# Validation predictions
pikc_pred <- read.csv('../../../Pikh1EngineeringProject/3_Full_Modeling/Output/Avr-PikC/esm2_t6_8M_UR50D_regression_predictions.csv')
pikf_pred <- read.csv('../../../Pikh1EngineeringProject/3_Full_Modeling/Output/Avr-PikF/esm2_t6_8M_UR50D_regression_predictions.csv')
pikc_pred <- pikc_pred %>% 
  mutate("Binding" = ifelse(enrichment > 0, "Positive", "Negative"))
pikf_pred <- pikf_pred %>% 
  mutate("Binding" = ifelse(enrichment > 0, "Positive", "Negative"))
predictions <- rbind(pikc_pred, pikf_pred)
# Scatter plots 
## Avr-PikC
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
ggsave('../Output/Figures/Fig1b.png', pikc_pred_plot, height = 56, width = 87, units="mm")
## Avr-PikF
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
ggsave('../Output/Figures/Fig1c.png', pikf_pred_plot, height = 56, width = 87, units="mm")

#################################################################

### Figure 2: MSA and clustering of sequences detected in the 3k RGP dataset
plot_fig2 <- function(a, b, c){
  # Dendogram
  myseqs <- read.alignment(paste0('../Output/Data/', a, '_HMA_Variants_Unique_Aligned_Formatted.fa'), format = 'fasta')
  myseqs$nam <- gsub('_\\(..\\)', '', myseqs$nam)
  myseqs$nam <- gsub('_\\(.\\)', '', myseqs$nam)
  mat <- dist.alignment(myseqs, matrix = 'identity', gap=TRUE)
  hc <- hclust(mat, 'average')
  dd <- as.dendrogram(hc) 
  dd2 <- rotate(dd, rev(b))
  ## Overall clusters with labels
  rep_dendo <- ggdendrogram(dd2, rotate = TRUE, size = 2)
  ## Save dendogram (flipped) without labels
  ddata <- dendro_data(dd2, type = "rectangle")
  deno <- ggplot(segment(ddata)) + 
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend), size=0.65) +
    coord_flip() +
    scale_y_reverse(expand = c(0.01, 0)) +
    scale_x_continuous(limits = c(-1.7, 14)) +
    theme_void() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank())
  ggsave(paste0('../Output/Figures/Fig2_', a, '_Dendo.png'), deno, height = 2.25, width = 1)
  # Sequence alignment
  msa <- ggmsa(paste0('../Output/Data/', a, '_HMA_Variants_Unique_Aligned_Formatted.fa'), 1, 82, seq_name=TRUE, char_width = 0.55, color = "Zappo_AA", border = "white", consensus_views = TRUE, use_dot = TRUE, ref = c) + 
    #geom_seqlogo(color = "Zappo_AA") +
    scale_x_continuous(breaks = seq(0, 82, by = 10), expand = c(0, 0)) +
    scale_y_discrete(labels=function(x) gsub("_", " ", x, fixed=TRUE)) +
    theme(text=element_text(size=5.25))
  ggsave(paste0('../Output/Figures/Fig2b_', a, '.png'), msa, height = 31.5, width = 158, units="mm", dpi=1000)
  return(msa)
}
pikp1_dendo <- plot_fig2('Pikp1', pikp_seq_order, "Pikh-1_(24)")
pikm1_dendo <- plot_fig2('Pikm1', pikm_seq_order, "Pikm-1_(14)")
