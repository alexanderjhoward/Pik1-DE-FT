# Script setup
library(tidyverse)
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
OUT_DIR <- "../Output/YSD"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
REFERENCE <- "Pikh-1"
SIGNAL_STAT <- "mean"

## Step 1: Read and format data
FILES <- list.files("../Data/YSD", pattern = "^per_cell\\.csv$", recursive = TRUE, full.names = TRUE)
names(FILES) <- basename(dirname(dirname(FILES)))
GENOTYPE_LEVELS <- c("pikh"="Pikh-1", "vel"="VK", "san"="SHZ-2")
LIGAND_LEVELS  <- c("a" = "Avr-PikA", "c" = "Avr-PikC")
cells <- read_csv(FILES, show_col_types = FALSE) %>%
  mutate(base = gsub("-S", "", base)) %>%
  mutate(session = folder) %>%
  separate(base, c("genotype", "ligand"), sep="_") %>%
  mutate(genotype = tolower(genotype),
         ligand = tolower(ligand)) %>%
  mutate(genotype = gsub("vk", "vel", genotype),
         genotype = gsub("szh2", "san", genotype)) %>%
  mutate(genotype = recode(genotype, !!!GENOTYPE_LEVELS),
         ligand = recode(ligand, !!!LIGAND_LEVELS)) %>% 
  mutate(genotype = factor(genotype, levels = GENOTYPE_LEVELS),
         ligand = factor(ligand, levels = LIGAND_LEVELS),
         session = factor(session))

## Step 2: Correct measurements against background (median out-of-cell signal)
stat_col <- function(ch) paste0(ch, "_", SIGNAL_STAT)
cells <- cells %>%
  mutate(expr_sig = .data[[stat_col("expr")]] - expr_bg,
         ligand_sig = .data[[stat_col("ligand")]] - ligand_bg) %>%
  filter(expr_sig > 0, ligand_sig > 0) # Keep anything above background staining

## Step 3: Take log2 ratio of filtered ligand/expression measurements
expr_floor <- cells %>%
  group_by(session, genotype) %>%
  summarise(expr_floor = quantile(expr_sig, 0.75, na.rm = TRUE), .groups = "drop") %>% #Only top 25% of receptor expressors per session within genotype are used
  mutate(merge_id = paste(session, genotype, sep="_")) %>%
  select(merge_id, expr_floor)
cells_pos <- cells %>%
  mutate(merge_id = paste(session, genotype, sep="_")) %>%
  left_join(expr_floor, by = "merge_id") %>%
  filter(expr_sig > expr_floor) %>% # Select cells by expression floor
  mutate(norm = log2(ligand_sig / expr_sig))

# Step 4: Normalize scores relative to Pikh-1
ref_level <- cells_pos %>%
  filter(genotype == REFERENCE) %>%
  group_by(session, ligand) %>%
  summarise(ref_norm = mean(norm), .groups = "drop")
cells_rel <- cells_pos %>%
  left_join(ref_level, by = c("session", "ligand")) %>%
  mutate(rel = norm - ref_norm)

# Save stats on normalized log2 ligand/receptor ratios
q_df <- cells_rel %>%
  group_by(ligand, genotype) %>%
  summarise(min = min(rel),
            q10 = quantile(rel, 0.10),
            q25 = quantile(rel, 0.25),
            med = median(rel),
            q75 = quantile(rel, 0.75), 
            q90 = quantile(rel, 0.90),
            max = max(rel),
            mean_rel = mean(rel),
            n_cells = n(), 
            .groups = "drop")
plot_df <- cells_rel %>%
  group_by(session, genotype, ligand) %>%
  summarise(rel = mean(rel), n_cells = n(), .groups = "drop")
write_csv(q_df, file.path(OUT_DIR, "ysd_quantiles.csv"))
write_csv(plot_df, file.path(OUT_DIR, "ysd_session_means.csv"))
