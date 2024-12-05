# Set up
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
library(tidyverse)

### Format sequence variants of the Pik-1 ligand binding domain

# Load in sample annotations
ids <- read.csv('../Source/accession_info.csv')
annots <- read.table('../Source/3K_list_sra_ids.txt', sep='\t', header=T) %>%
  mutate(Country_Origin_updated = gsub(" ", "'", Country_Origin_updated)) %>%
  mutate(Country_Origin_updated = gsub("Korea,_Republic_of", "Republic_of_Korea", Country_Origin_updated)) %>%
  mutate(Country_Origin_updated = gsub("_no_info", "Not_Available", Country_Origin_updated))
sample_info <- merge(annots, ids, by.x="SRA.Accession", by.y="SRA")
rm(ids, annots)

# Load in and format exons
exon1 <- read.csv(paste0('../Output/N22/N22_E1_pep.csv'))
exon1 <- separate(exon1, id, into=c('BIOSAMPLE', "Exon1"), sep="_")
exon2 <- read.csv(paste0('../Output/N22/N22_E2_pep.csv'))
exon2 <- separate(exon2, id, into=c('BIOSAMPLE', "Exon2"), sep="_")

# Merge exons into contiguous sequence
pik <- merge(sample_info, exon1, by='BIOSAMPLE')
pik <- merge(pik, exon2, by='BIOSAMPLE')
rm(exon1, exon2)

# Remove heterozygous reads identified manually in a genome browser
het_reads <- c('GLKQKIVFKIPMEGNNCRSKAMALVASTGGVDSVALVGDLRDKIEVVGDGIDSINLVSALRKKVGPAMFLEVSQAKED',
               'GLKQTIVIKVAMEGNNCRSKAMDLVKSTGGVYSVSLAGDLRDKIEVVGYGIDPIKLISALRKKVGHAELLQVSQAKKD')
pik <- pik %>%
  mutate(SEQ = paste(seq.x, seq.y, sep="")) %>%
  filter(!(SEQ %in% het_reads)) %>%
  select(c("BIOSAMPLE", "SRA.Accession", "X3K_DNA_IRIS_UNIQUE_ID", "Genetic_Stock_varname", "Country_Origin_updated", "SEQ"))
colnames(pik) <- c("BIOSAMPLE", "SRA", "IRIS", "STOCK", "COUNTRY", "SEQ")
rm(het_reads)

# Back-fill the C-terminal deletion found in 2 samples
pik <- pik %>% 
  mutate(SEQ = ifelse(SEQ == "GLKQKIVIKVAMEGNNCRSKAMALVASTGGVDSVALVGDLRDKIVVVGYGIDPIKLISALRKKVGDAELLQVSD", paste0(SEQ, "VKET"), SEQ))

# Save sequence results for all searched varieties
write.csv(pik, file = paste0('../Output/Pikp_HMA_Variants.csv'),  row.names = FALSE)

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
rm(pik_freq)

# Format names and save all unique variants
pik_uniq <- pik_uniq %>%
  mutate(Variety = gsub("Iet_14720", "Pikh-1", Variety)) %>%
  mutate(Variety = gsub("Simul_Khuri", "Pikp-1", Variety))
write.csv(pik_uniq, file = '../Output/Pikp_HMA_Variants_Unique.csv',  row.names = FALSE)

# Save FASTA of unique variants in specified order (optional ordering, just for aesthetics)
seq_order <- c("Pikh-1", "Pikp-1", "Bangla", "Larha_Mugad", "Chundi", 
               "Periya_Vellai", "Vellai_Kolomban", "Taipei_Woo_Co", "Adukkan", 
               "Xintuan_Hei_Gu", "Zimangfeie", "Banjiemang2", "Sanhuangzhan_No_2")
fasta <- pik_uniq %>%
  mutate(Variety = factor(Variety, levels = seq_order)) %>%
  arrange(Variety) %>%
  mutate(Number = sprintf("%02d", 1:nrow(pik_uniq))) %>%
  mutate(Variety = paste0('>', Number, ";", Variety, "_(", Freq, ")")) %>%
  select(Variety, SEQ)
fasta <- do.call(rbind, lapply(seq(nrow(fasta)), function(i) t(fasta[i, ])))
write.table(fasta, file = '../Output/Pikp_HMA_Variants_Unique.fa',  row.names = FALSE, col.names = FALSE, quote = FALSE)

#######

### Format sequence variants of the full-length Pik-1 and Pik-2 sequences containing SHZ-2 or VK ligand binding domains

# Load in and format Pik-1 and Pik-2 exons
p1_exon1 <- read.csv('../Output/Vars/Vars_P1E1_pep.csv')
p1_exon1 <- separate(p1_exon1, id, into=c('BIOSAMPLE', "Exon1"), sep="_")
p1_exon2 <- read.csv('../Output/Vars/Vars_P1E2_pep.csv')
p1_exon2 <- separate(p1_exon2, id, into=c('BIOSAMPLE', "Exon2"), sep="_")
p1_exon3 <- read.csv('../Output/Vars/Vars_P1E3_pep.csv')
p1_exon3 <- separate(p1_exon3, id, into=c('BIOSAMPLE', "Exon3"), sep="_")
p2_exon1 <- read.csv('../Output/Vars/Vars_P2E1_pep.csv')
p2_exon1 <- separate(p2_exon1, id, into=c('BIOSAMPLE', "Exon1"), sep="_")
p2_exon2 <- read.csv('../Output/Vars/Vars_P2E2_pep.csv')
p2_exon2 <- separate(p2_exon2, id, into=c('BIOSAMPLE', "Exon2"), sep="_")

# Merge exons into contiguous sequences
p1 <- merge(sample_info, p1_exon1, by='BIOSAMPLE')
p1 <- merge(p1, p1_exon2, by='BIOSAMPLE')
p1 <- merge(p1, p1_exon3, by='BIOSAMPLE')
p1$SEQ <- paste(p1$seq.x, p1$seq.y, p1$seq, sep="")
p2 <- merge(sample_info, p2_exon1, by='BIOSAMPLE')
p2 <- merge(p2, p2_exon2, by='BIOSAMPLE')
p2$SEQ <- paste(p2$seq.x, p2$seq.y, sep="")
rm(p1_exon1, p1_exon2, p1_exon3, p2_exon1, p2_exon2)

# Format data
p1 <- p1 %>%
  select(c("BIOSAMPLE", "SRA.Accession", "X3K_DNA_IRIS_UNIQUE_ID", "Genetic_Stock_varname", "Country_Origin_updated", "SEQ"))
colnames(p1) <- c("BIOSAMPLE", "SRA", "IRIS", "STOCK", "COUNTRY", "SEQ")
p2 <- p2 %>%
  select(c("BIOSAMPLE", "SRA.Accession", "X3K_DNA_IRIS_UNIQUE_ID", "Genetic_Stock_varname", "Country_Origin_updated", "SEQ"))
colnames(p2) <- c("BIOSAMPLE", "SRA", "IRIS", "STOCK", "COUNTRY", "SEQ")

# Save sequence results for all searched varieties
write.csv(p1, file = '../Output/Vars_Pik1_Sequences.csv',  row.names = FALSE)
write.csv(p2, file = '../Output/Vars_Pik2_Sequences.csv',  row.names = FALSE)

# Identify frequency of variants
p1_freq <- as.data.frame(table(p1$SEQ))
p2_freq <- as.data.frame(table(p2$SEQ))

# Filter for unique variants and save
p1_uniq <- p1 %>% distinct(SEQ, .keep_all = TRUE)
p2_uniq <- p2 %>% distinct(SEQ, .keep_all = TRUE)
p1_uniq <- merge(p1_uniq, p1_freq, by.x="SEQ", by.y="Var1") %>% 
  separate(STOCK, c("Variety", "ID"), sep="::") %>%
  mutate(Variety = gsub("(?<=\\p{L})(\\p{L}+)", "\\L\\1", Variety, perl = TRUE))
p2_uniq <- merge(p2_uniq, p2_freq, by.x="SEQ", by.y="Var1") %>% 
  separate(STOCK, c("Variety", "ID"), sep="::") %>%
  mutate(Variety = gsub("(?<=\\p{L})(\\p{L}+)", "\\L\\1", Variety, perl = TRUE))
write.csv(p1_uniq, file = '../Output/Vars_Pik1_Sequences_Unique.csv',  row.names = FALSE)
write.csv(p2_uniq, file = '../Output/Vars_Pik2_Sequences_Unique.csv',  row.names = FALSE)
