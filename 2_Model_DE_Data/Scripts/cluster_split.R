# Set up
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
library(tidyverse)
library(plotly)

# Use tSNE dimensionality reduction on Avr-PikC and Avr-PikF sequence embeddings
library(Rtsne)
calc_tsne <- function(a){
  emb <- a
  matrix <- as.matrix(emb[,3:322])
  set.seed(42)
  tsne <- Rtsne(matrix)
  tsne_df <- data.frame(hma = emb$hma,
                        enrichment = emb$enrichment,
                        x = tsne$Y[,1], 
                        y = tsne$Y[,2])
  return(tsne_df)
}
pikc_emb <- rbind(read.csv(file='../Output/AvrPikC/train_embeddings.csv'),
                  read.csv(file='../Output/AvrPikC/val_embeddings.csv'))
pikc_tsne_df <- calc_tsne(pikc_emb)
pikf_emb <- rbind(read.csv(file='../Output/AvrPikF/train_embeddings.csv'),
                  read.csv(file='../Output/AvrPikF/val_embeddings.csv'))
pikf_tsne_df <- calc_tsne(pikf_emb)

# Plot tSNE results
pikc_tsne_plot <- ggplotly(ggplot(pikc_tsne_df) + 
                             geom_point(aes(x=x,y=y,color=enrichment)) +
                             scale_colour_gradient2()+
                             theme_bw())
pikf_tsne_plot <- ggplotly(ggplot(pikf_tsne_df) + 
                             geom_point(aes(x=x,y=y,color=enrichment)) +
                             scale_colour_gradient2() +
                             theme_bw())

# Train/validation split Avr-PikC and Avr-PikF data by cluster (trying to maintain score distributions and 0.95/0.05 ratio)
## Avr-PikC
pikc_cluster1 <- pikc_tsne_df %>%
  filter((x > 40) & (y < -10))
pikc_cluster2 <- pikc_tsne_df %>%
  filter((x > -5) & (y > 36))
pikc_cluster3 <- pikc_tsne_df %>%
  filter(x < -45)
pikc_cluster <- rbind(pikc_cluster1, pikc_cluster2, pikc_cluster3)
rm(pikc_cluster1, pikc_cluster2, pikc_cluster3)
pikc_emb_train <- pikc_emb %>%
  filter(!(hma %in% pikc_cluster$hma))
pikc_emb_val <- pikc_emb %>%
  filter(hma %in% pikc_cluster$hma)
ggplot(pikc_emb_val, aes(enrichment)) + geom_histogram()
ggplot(pikc_emb_train, aes(enrichment)) + geom_histogram()
write.csv(pikc_emb_train, '../Output/AvrPikC/train_embeddings_cluster_split.csv', row.names = F)
write.csv(pikc_emb_val, '../Output/AvrPikC/val_embeddings_cluster_split.csv', row.names = F)

## Avr-PikF
pikf_cluster1 <- pikf_tsne_df %>%
  filter(y > 45)
pikf_cluster2 <- pikf_tsne_df %>%
  filter(x < -45)
pikf_cluster3 <- pikf_tsne_df %>%
  filter((x > 0) & (x < 10) & (y < -35))
pikf_cluster <- rbind(pikf_cluster1, pikf_cluster2, pikf_cluster3)
rm(pikf_cluster1, pikf_cluster2, pikf_cluster3)
pikf_emb_train <- pikf_emb %>%
  filter(!(hma %in% pikf_cluster$hma))
pikf_emb_val <- pikf_emb %>%
  filter(hma %in% pikf_cluster$hma)
ggplot(pikf_emb_val, aes(enrichment)) + geom_histogram()
ggplot(pikf_emb_train, aes(enrichment)) + geom_histogram()
write.csv(pikf_emb_train, '../Output/AvrPikF/train_embeddings_cluster_split.csv', row.names = F)
write.csv(pikf_emb_val, '../Output/AvrPikF/val_embeddings_cluster_split.csv', row.names = F)
