library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
library(tidyverse)
library(readxl)

# Function for cleaning directed evolution sequencing data
clean_facs <- function(x, y){
  
  ## Read in raw sequencing data
  raw_reads <- read_excel(x)
  message(paste('Raw read row count:', nrow(raw_reads)))
  
  ## Remove too short or too long protein sequences
  raw_reads$seqlen <- apply(raw_reads,2,nchar)[,7]
  summary(raw_reads$seqlen)
  message(paste('Truncated sequences removed:', nrow(raw_reads %>% filter(seqlen < 99))))
  message(paste('Extended sequences removed:', nrow(raw_reads %>% filter(seqlen > 99)) ))
  filtered_1 <- raw_reads %>% filter(seqlen == 99)
  message(paste('99 aa filter row count:', nrow(filtered_1)))

  ## Find and remove sequences containing "X" amino acids
  message(paste('Ambiguous sequences removed:', nrow(filtered_1 %>% filter(grepl("X", TargetAA)))))
  filtered_2 <- filtered_1 %>% filter(!(grepl("X", TargetAA)))
  message(paste('X aa filter row count:', nrow(filtered_2)))
  
  ## Separate out linker sequences from HMA domain
  filtered_2$nlink <- substr(filtered_2$TargetAA,1,19)
  filtered_2$hma <- substr(filtered_2$TargetAA,20,97)
  filtered_2$clink <- substr(filtered_2$TargetAA,98,99)
  
  ## Combine identical HMA domains and sum up counts
  final <- filtered_2 %>%
    group_by(hma) %>%
    summarize(reads = sum(Reads)) %>%
  return(final)
  
}

# Clean sequencing data
lib <- clean_facs('../Source/20231120_library1_raw.xlsx')
sort1c <- clean_facs('../Source/20231120_sort1_C_raw.xlsx')
sort1f <- clean_facs('../Source/20231120_sort1_F_raw.xlsx')

# Function for calculating enrichment score based on read count
combine_seqs <- function(a, b, n){
  seqs_of_int <- c("GLKQKIVIKVAMEGNNCRSKAMALVASTGGVDSVALVGDLRDKIEVVGYGIDPIKLISALRKKVGDAELLQVSQAKKD", # WT
                   "GLKRIIVIKVAREGNNCRSKAMALVASTGGVDSVALVGDLRDKIEVVGYGIDPIKLISALRKKVGDAELLQVSQAKKD", # 5_1, 6_1
                   "GLKRKIVIKVAVEGNNCRSKAMALVASTGGVDSVALVGDLRDKIEVVGCGIDPIKLISALRKKVGGAVLLQISQAKKD", # 5_2, 6_3
                   "GLKQKIVIKVAMEGNNCRSKAMALVASTGGVDSVALVGDLRDKIEVVGYGIDPTKLISALRKKVGGAELLQVSRVKKG", # FCFC 11 (new name FCFC1=B4.1)
                   "GLKQKIVIKVAMEGNNCRPKAMALVASTGGVESAALVGDLRDKIEVVGHGIDPIKLISALRKKVGDAELSQVSLVKKD", # FCFC 6 (new name FCFC2)
                   "GLKQKIVIKVAMEGNNCRSKAMALVASTGGVDSVALVGDLWDKIEVVGYGIDLIKLISALRKKVGDAELLQVSQVKKG", # CFCF 2 (new name CFCF1)
                   "GLKQKIVIKVAKEGNDCRSRAMALVASTGGVDSVALVGDLRDKIEVVGYGIDPIKLISALRKKVGDAELLQVSRVKKD", # CFCF 7 (new name CFCF2=A4.1)
                   "GLKRIIVIKVARGGNNCRSKAMALVASTGGVDSVALVGDLRDKIEVVGYGIDPIKLISALRKKVGDAELLQVSQAKKG") # 7.12 (new name 7.1)
  combo <- merge(a, b, by="hma", all=T)
  combo$reads.x[is.na(combo$reads.x)] <- 1 # Sequences missing in the starting library but present in the ending library must be made present in the starting library (+1 read)
  combo[is.na(combo)] <- 0 # Fill in any remaining NA values with 0 reads
  combo <- combo %>% 
    mutate(norm_reads.x = (reads.x)/sum(reads.x)) %>% # Proportionally normalize reads within starting library
    mutate(norm_reads.y = (reads.y+1)/sum(reads.y)) %>% # Proportionally normalize reads within post-sort library, adding 1 to the numerator to avoid zeros
    mutate(enrichment = log10(norm_reads.y/norm_reads.x)) %>% # Log normalize data
    mutate(present = ifelse(enrichment <= 0, 0, 1)) %>% 
    filter(!(reads.x <= n & reads.y <= n)) %>% # After calculating scores, only keep sequences that have more than n reads in the starting or ending library
    filter(!(hma %in% seqs_of_int)) %>%
    mutate(bin = ntile(enrichment, round(n()/100, digits = 0)))
  return(combo)
}

# Calculate enrichment scores for each ligand
pikc <- combine_seqs(lib, sort1c, 9) # Avr-PikC binding at 1 uM (only including sequences with at least 10 reads in either library)
pikf <- combine_seqs(lib, sort1f, 9) # Avr-PikF binding at 1 uM (only including sequences with at least 10 reads in either library)

# Save enrichment data
write.csv(pikc, file='../Output/Data/AvrPikC_Enrichment.csv', row.names = F)
write.csv(pikf, file='../Output/Data/AvrPikF_Enrichment.csv', row.names = F)

# Plot distribution of Avr-PikC and Avr-PikF enrichment scores
pikc <- pikc %>% mutate(Ligand='Avr-PikC')
pikf <- pikf %>% mutate(Ligand='Avr-PikF')
ligands <- rbind(pikc, pikf)
ligands_plot <- ggplot(ligands, aes(x=enrichment, fill=Ligand)) +
  geom_histogram(color='black') +
  geom_vline(xintercept=0, color='red', linetype='dashed') +
  scale_fill_manual(values=c('#E1AA7D', 'lightblue')) +
  facet_wrap(~Ligand) +
  xlab('Enrichment Score') +
  ylab('Read Count') +
  theme_minimal()
ggsave('../Output/Plots/enrichment_distribution_avrpikc_avrpikf.jpg', ligands_plot, height = 2, width=6)
