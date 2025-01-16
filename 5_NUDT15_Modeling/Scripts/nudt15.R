# Set up
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
library(tidyverse)
library(seqinr)

# Load and save  wt information
wt <- read.csv(file='../Source/NUDT15.csv') %>%
  filter(Info=="WT")
wtseq <- wt$Protein[1]
rm(wt)

# List of variants with known clinical phenotypes collected from literature search (found in Source directory)
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
genomad <- read.csv('../Source/gnomAD_v4.1.0_ENSG00000136159_2025_01_13_14_21_37.csv') %>%
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
df <- read.csv(file='../Source/NUDT15.csv') %>%
  drop_na(Final.NUDT15.activity.Score) %>%
  mutate(Score=Final.NUDT15.activity.Score) %>%
  select(Protein, Score, Final_classification) %>% # After calculating scores, only keep sequences that have more than n reads in the starting or ending library
  mutate(Score = log(Score + 0.55)) %>% # Minimum cut-off for wt-like activity set by authors was 0.45, so add 0.55 to make all wt-like activity >= 1
  filter(!(Protein %in% variants$seq)) %>%
  mutate(bin = ntile(Score, round(n()/100, digits = 0)))
write.csv(df, file='../Output/Data/NUDT15_cleaned.csv', row.names = F)

# Remove any genomAD variants that have been already tested in this assay
genomad_vars <- genomad_vars %>%
  filter(!(seq %in% df$Protein))
full_vars <- rbind(variants, genomad_vars)
rm(variants, genomad_vars)
write.csv(full_vars, file='../Output/Data/NUDT15_variants.csv', row.names = F)
