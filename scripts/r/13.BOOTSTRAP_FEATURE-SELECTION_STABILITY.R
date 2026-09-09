# BOOTSTRAP FEATURE-SELECTION STABILITY (TRAIN samples only)

library(limma)
library(AnnotationDbi)      
library(hgu133plus2.db)     
library(ggplot2)            


expr_train <- readRDS("expr_train.rds")
pheno_clean <- readRDS("pheno_clean.rds")

common_genes <- readRDS("common_genes.rds")


stage_order <- c(
  "Normal Colon",
  "Adenoma (Polyp)",
  "Primary Colorectal Cancer",
  "Metastasis (Liver/Lung)"
)


fitStep_fast <- function(x) {
  n <- length(x)
  cs  <- cumsum(x)
  css <- cumsum(x^2)
  total   <- cs[n]
  totalsq <- css[n]
  
  i  <- 2:(n - 2)
  n1 <- i
  n2 <- n - i
  s1 <- cs[i];  s2 <- total - s1
  ss1 <- css[i]; ss2 <- totalsq - ss1
  sse <- (ss1 - s1^2 / n1) + (ss2 - s2^2 / n2)
  
  best <- which.min(sse)
  best_index <- i[best]
  best_sse   <- sse[best]
  sstot <- totalsq - total^2 / n
  statistic <- ((sstot - best_sse) / 3) / (best_sse / (n - 4))
  
  list(index = best_index, statistic = statistic)
}




expr_probe_clean <- expr_train
if (max(expr_probe_clean, na.rm = TRUE) > 50) {
  expr_probe_clean <- log2(expr_probe_clean + 1)
}
keep_genes_probe <- rowMeans(expr_probe_clean) > 5
expr_probe_clean <- expr_probe_clean[keep_genes_probe, ]

cat("Probe-level train matrix rebuilt:", nrow(expr_probe_clean), "probes x",
    ncol(expr_probe_clean), "samples\n")


set.seed(123)
n_boot <- 50   

samples_train <- rownames(pheno_clean)
stage_train   <- pheno_clean$Stage


probe_symbol_map <- mapIds(
  hgu133plus2.db, keys = rownames(expr_probe_clean),
  column = "SYMBOL", keytype = "PROBEID", multiVals = "first"
)

boot_gene_hits <- character(0)   # accumulates gene symbols selected each bootstrap

for (b in seq_len(n_boot)) {
  
  # Stratified bootstrap: resample WITH replacement, separately per stage
  boot_idx <- unlist(lapply(levels(stage_train), function(s) {
    idx <- which(stage_train == s)
    sample(idx, size = length(idx), replace = TRUE)
  }))
  
  boot_ids   <- samples_train[boot_idx]
  boot_cols  <- make.unique(boot_ids)     # unique column names for repeats
  
  expr_b  <- expr_probe_clean[, boot_ids]  
  colnames(expr_b) <- boot_cols
  pheno_b <- pheno_clean[boot_ids, , drop = FALSE]
  rownames(pheno_b) <- boot_cols
  pheno_b$Stage <- factor(pheno_b$Stage, levels = stage_order)
  
  # LIMMA on this bootstrap resample 
  
  design_b <- model.matrix(~0 + Stage, data = pheno_b)
  colnames(design_b) <- make.names(levels(pheno_b$Stage))
  contrast_b <- makeContrasts(
    Adenoma_vs_Normal    = Adenoma..Polyp. - Normal.Colon,
    Cancer_vs_Adenoma    = Primary.Colorectal.Cancer - Adenoma..Polyp.,
    Metastasis_vs_Cancer = Metastasis..Liver.Lung. - Primary.Colorectal.Cancer,
    levels = design_b
  )
  fit_b  <- lmFit(expr_b, design_b)
  fit2_b <- eBayes(contrasts.fit(fit_b, contrast_b))
  
  deg_b <- c()
  for (comp in colnames(contrast_b)) {
    res_b <- topTable(fit2_b, coef = comp, number = Inf, adjust.method = "BH")
    deg_b <- c(deg_b, rownames(res_b)[res_b$adj.P.Val < 0.05 & abs(res_b$logFC) > 1])
  }
  deg_b <- unique(deg_b)
  if (length(deg_b) < 10) next
  
  #  monotonic filter 
  stage_mean_b <- sapply(levels(pheno_b$Stage), function(s) {
    ids <- rownames(pheno_b)[pheno_b$Stage == s]
    rowMeans(expr_b[deg_b, ids, drop = FALSE])
  })
  
  is_monotonic <- function(x) all(diff(x) >= 0) || all(diff(x) <= 0)
  genes_mono_b <- deg_b[apply(stage_mean_b, 1, is_monotonic)]
  if (length(genes_mono_b) < 10) next
  
  # StepMiner (fast vectorized version, ordered by stage) 
  
  ordered_b <- rownames(pheno_b)[order(pheno_b$Stage)]
  expr_ord_b <- expr_b[genes_mono_b, ordered_b, drop = FALSE]
  
  step_b <- data.frame()
  for (g in genes_mono_b) {
    r <- fitStep_fast(as.numeric(expr_ord_b[g, ]))
    step_b <- rbind(step_b, data.frame(ProbeID = g, StepIndex = r$index, F_Stat = r$statistic))
  }
  step_b <- step_b[step_b$F_Stat > 10, ]
  if (nrow(step_b) == 0) next
  step_b$SYMBOL <- probe_symbol_map[step_b$ProbeID]
  step_b <- step_b[!is.na(step_b$SYMBOL), ]
  
  # boundaries + transition assignment for THIS bootstrap's ordering 
  
  stage_levels_b <- levels(pheno_b$Stage)
  smap_b <- data.frame(SampleOrder = seq_along(ordered_b), SampleID = ordered_b,
                       Stage = pheno_b[ordered_b, "Stage"])
  bounds_b <- data.frame()
  for (i in 1:(length(stage_levels_b) - 1)) {
    cur <- stage_levels_b[i]; nxt <- stage_levels_b[i + 1]
    last_idx <- max(smap_b$SampleOrder[smap_b$Stage == cur])
    bounds_b <- rbind(bounds_b, data.frame(Boundary_Index = last_idx,
                                           Transition = paste(cur, "→", nxt)))
  }
  assign_transition_b <- function(idx) {
    for (i in 1:nrow(bounds_b)) if (idx <= bounds_b$Boundary_Index[i]) return(bounds_b$Transition[i])
    bounds_b$Transition[nrow(bounds_b)]
  }
  step_b$Transition <- sapply(step_b$StepIndex, assign_transition_b)
  
  # top-25-per-transition (same rule as our main pipeline)
  
  top25_b <- unique(unlist(lapply(unique(step_b$Transition), function(tr) {
    sub <- step_b[step_b$Transition == tr, ]
    sub <- sub[!duplicated(sub$SYMBOL), ]
    sub <- sub[order(-sub$F_Stat), ]
    head(sub$SYMBOL, 25)
  })))
  
  boot_gene_hits <- c(boot_gene_hits, top25_b)
  
  if (b %% 10 == 0) cat("  bootstrap", b, "/", n_boot, "done\n")
}


# SUMMARIZE STABILITY


freq_table <- table(boot_gene_hits)
freq_df <- data.frame(SYMBOL = names(freq_table), Times_Selected = as.integer(freq_table))
freq_df$Selection_Rate <- freq_df$Times_Selected / n_boot
freq_df <- freq_df[order(-freq_df$Selection_Rate), ]

write.csv(freq_df, "Bootstrap_Gene_Stability_All.csv", row.names = FALSE)

# focus on OUR 72 final-panel genes specifically

final_panel_stability <- freq_df[freq_df$SYMBOL %in% common_genes, ]
missing_from_boot <- setdiff(common_genes, freq_df$SYMBOL)
if (length(missing_from_boot) > 0) {
  final_panel_stability <- rbind(final_panel_stability,
                                 data.frame(SYMBOL = missing_from_boot, Times_Selected = 0, Selection_Rate = 0))
}
final_panel_stability <- final_panel_stability[order(-final_panel_stability$Selection_Rate), ]
write.csv(final_panel_stability, "Bootstrap_Gene_Stability_FinalPanel.csv", row.names = FALSE)


cat("\nBOOTSTRAP STABILITY SUMMARY (", n_boot, "resamples):\n")
cat("Mean selection rate across our 72-gene panel:",
    round(mean(final_panel_stability$Selection_Rate) * 100, 1), "%\n")
cat("Genes selected in >=80% of bootstraps:",
    sum(final_panel_stability$Selection_Rate >= 0.8), "/", nrow(final_panel_stability), "\n")
cat("Genes selected in <20% of bootstraps (unstable):",
    sum(final_panel_stability$Selection_Rate < 0.2), "/", nrow(final_panel_stability), "\n")
print(head(final_panel_stability, 15))




#Figure

top20_stability <- head(final_panel_stability, 20)
top20_stability$SYMBOL <- factor(top20_stability$SYMBOL, levels = rev(top20_stability$SYMBOL))
p_s3 <- ggplot(top20_stability, aes(x = SYMBOL, y = Selection_Rate)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(title = "Bootstrap Gene-Selection Stability (Top 20 of 72-gene panel)",
       x = "", y = "Selection Rate (50 bootstraps)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
print(p_s3)
ggsave("FigureS3_Bootstrap_Gene_Stability.png", p_s3, width = 8, height = 6, dpi = 300)

#Figure End
