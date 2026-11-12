# ---------------------------------------------------------------------------
# if (!requireNamespace("BiocManager")) install.packages("BiocManager")
#   BiocManager::install(c("flowCore", "flowStats"))
#   install.packages(c("tidyverse", "ggridges"))
# ---------------------------------------------------------------------------
library(flowCore)
library(flowStats)
library(tidyverse)
library(ggridges)
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

# Step 1: Set up script and load FACS data
FCS_DIR <- "../Source/FACS"
OUT_DIR <- "../Output/FACS"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CH_EXPR   <- "Alexa Fluor 488-A"
CH_LIGAND <- "PE-A"
CH_FSC_H  <- "FSC-H"
CH_SSC_H  <- "SSC-H"
CH_SSC_W  <- "SSC-W"
GENOTYPE_LEVELS <- c("Pikh-1", "VK", "SHZ-2")

meta <- tribble(~file, ~genotype, ~ligand, ~replicate, ~role,
                "C__H.fcs", "Pikh-1", "Avr-PikC", 1L, "sample",
                "C__V.fcs", "VK", "Avr-PikC", 1L, "sample",
                "C__S.fcs", "SHZ-2", "Avr-PikC", 1L, "sample",
                "A__H.fcs", "Pikh-1", "Avr-PikA", 1L, "sample",
                "A__V.fcs", "VK", "Avr-PikA", 1L, "sample",
                "A__S.fcs", "SHZ-2", "Avr-PikA", 1L, "sample",
                "Controls_H__D16_UNS.fcs", NA, NA, 1L, "unstained",
                "Controls_H__D16_488Only.fcs", NA, NA, 1L, "comp_expr",
                "Controls_H__D16_PE_Only.fcs", NA, NA, 1L, "comp_ligand",)

fs <- read.flowSet(files = file.path(FCS_DIR, meta$file),
                   transformation = FALSE,
                   truncate_max_range = FALSE)
sampleNames(fs) <- meta$file

# Step 1: Compensate for spillover between between FITC and PE stains
get_ff <- function(role) {
  f <- meta$file[meta$role == role & !is.na(meta$role)]
  if (length(f) == 0) return(NULL)
  fs[[f[1]]]
}
ctrl_expr <- get_ff("comp_expr")
ctrl_ligand <- get_ff("comp_ligand")
ctrl_unst <- get_ff("unstained")
af_expr <- median(exprs(ctrl_unst)[, CH_EXPR])
af_ligand <- median(exprs(ctrl_unst)[, CH_LIGAND])

bleed <- function(ff, from_ch, into_ch, af_from, af_into) {
  e <- exprs(ff)
  x <- e[, from_ch] - af_from # From emission channel (median normalized)
  y <- e[, into_ch] - af_into # Into measured channel (median normalized)
  pos <- x > quantile(x, 0.75, na.rm = TRUE)   # Select top 75% stained population
  median(y[pos], na.rm = TRUE) / median(x[pos], na.rm = TRUE)
}
s_e2l <- bleed(ctrl_expr, CH_EXPR, CH_LIGAND, af_expr, af_ligand)
s_l2e <- bleed(ctrl_ligand, CH_LIGAND, CH_EXPR, af_ligand, af_expr)
message(sprintf("Spillover %s -> %s: %.4f", CH_EXPR,   CH_LIGAND, s_e2l))
message(sprintf("Spillover %s -> %s: %.4f", CH_LIGAND, CH_EXPR,   s_l2e))
spill <- matrix(c(1, s_e2l,s_l2e, 1),
                nrow = 2, byrow = TRUE,
                dimnames = list(c(CH_EXPR, CH_LIGAND),
                                c(CH_EXPR, CH_LIGAND)))
write.csv(spill, file.path(OUT_DIR, "spillover_matrix.csv"))
fs <- compensate(fs, spill)

# Step 2: Set gates to select for single yeast cells
gate_one <- function(ff) {
  e <- exprs(ff)
  if (CH_FSC_H %in% colnames(ff)) {
    r <- e[, CH_FSC_H] / e[, CH_SSC_H]
    keep <- r > quantile(r, 0.05, na.rm = TRUE) &
    r < quantile(r, 0.95, na.rm = TRUE)
    ff <- ff[keep, ]
    }
  cells <- norm2Filter(CH_SSC_H, CH_SSC_W, scale.factor = 2.5)
  flowCore::Subset(ff, cells)
}
fs_gated <- fsApply(fs, gate_one)

# Step 3: Map FACS data with adjusted fluorescence, filter out poorly stained
#         cells, and take log2 of ratio between ligand and receptor signal
samples <- meta %>% dplyr::filter(role == "sample")
cells <- map_dfr(seq_len(nrow(samples)), function(i) {
  row <- samples[i, ]
  e <- exprs(fs_gated[[row$file]])
  tibble(
    file      = row$file,
    genotype  = row$genotype,
    ligand_id = row$ligand,
    replicate = row$replicate,
    expr      = e[, CH_EXPR] - af_expr,
    ligand    = e[, CH_LIGAND] - af_ligand
  )
}) %>%
  mutate(genotype = factor(genotype, levels = GENOTYPE_LEVELS))

trough_cut <- function(x) {
  x <- log10(pmax(x, 1))
  d <- density(x, n = 1024)
  lo <- which(d$x > quantile(x, 0.05) & d$x < quantile(x, 0.95))
  10^d$x[lo[which.min(d$y[lo])]]
}
floor_expr <- trough_cut(cells$expr)
message(sprintf("Expression-positive cutoff: %.1f", floor_expr))
pdf(file.path(OUT_DIR, "qc_expression_cutoff.pdf"), width = 5, height = 3.5)
plot(density(log10(pmax(cells$expr, 1))), main = "Expression, pooled samples",
     xlab = "log10 receptor expression")
abline(v = log10(floor_expr), col = "red", lty = 2)
dev.off()

cells_pos <- cells %>%
  dplyr::filter(expr > floor_expr, ligand > 0) %>%
  mutate(norm = log2(ligand / expr))
write_csv(cells_pos %>% count(file, genotype, replicate, name = "n_cells"),
          file.path(OUT_DIR, "cells_per_sample.csv"))
write_csv(cells_pos %>% select(genotype, ligand_id, norm),
          file.path(OUT_DIR, "facs_normalized_cells.csv"))
