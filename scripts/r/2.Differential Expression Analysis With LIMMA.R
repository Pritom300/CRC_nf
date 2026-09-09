library(limma)
library(GEOquery)
args <- commandArgs(trailingOnly = TRUE)
gse_file    <- args[1]
meta_file   <- args[2]
output_dir  <- args[3]


gse <- getGEO(filename = gse_file, getGPL = FALSE)
expr <- exprs(gse)


# Load metadata

pheno <- read.csv(meta_file)


rownames(pheno) <- pheno$geo_accession


# Match samples correctly

common_samples <- intersect(colnames(expr), rownames(pheno))

expr <- expr[, common_samples]
pheno <- pheno[common_samples, ]

cat("Matched samples:", length(common_samples), "\n")



# TRAIN/TEST SPLIT 

set.seed(123)

train_idx <- sample(seq_len(nrow(pheno)), size = floor(0.7 * nrow(pheno)))
train_samples <- rownames(pheno)[train_idx]
test_samples  <- rownames(pheno)[-train_idx]

pheno_train <- pheno[train_samples, ]
pheno_test  <- pheno[test_samples, ]


pheno_train <- pheno_train[pheno_train$Stage != "Unclassified/Other", ]
pheno_test  <- pheno_test[pheno_test$Stage != "Unclassified/Other", ]

expr_train <- expr[, rownames(pheno_train)]
expr_test  <- expr[, rownames(pheno_test)]

cat("Train samples:", ncol(expr_train), "| Test samples:", ncol(expr_test), "\n")




expr_clean <- expr_train
pheno_clean <- pheno_train




# Log2 check

if (max(expr_clean, na.rm = TRUE) > 50) {
  cat("Applying log2 transformation...\n")
  expr_clean <- log2(expr_clean + 1)
} else {
  cat("Already log2 scale\n")
}


# Filter genes

keep_genes <- rowMeans(expr_clean) > 5
expr_clean <- expr_clean[keep_genes, ]




# Design matrix

pheno_clean$Stage <- factor(pheno_clean$Stage)

design <- model.matrix(~0 + Stage, data = pheno_clean)
colnames(design) <- make.names(levels(pheno_clean$Stage))

print(colnames(design))


# Contrast Matrix 


colnames(design)
# check names manually once

contrast.matrix <- makeContrasts(
  Adenoma_vs_Normal = Adenoma..Polyp. - Normal.Colon,
  Cancer_vs_Adenoma = Primary.Colorectal.Cancer - Adenoma..Polyp.,
  Metastasis_vs_Cancer = Metastasis..Liver.Lung. - Primary.Colorectal.Cancer,
  levels = design
)


# LIMMA

fit <- lmFit(expr_clean, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


# DEG extraction
dir.create(file.path(output_dir, "LIMMA_Results"), showWarnings = FALSE)

for (comp in colnames(contrast.matrix)) {
  
  res <- topTable(fit2, coef = comp, number = Inf, adjust.method = "BH")
  
  res$ProbeID <- rownames(res)
  res$Comparison <- comp
  
  sig <- res[res$adj.P.Val < 0.05 & abs(res$logFC) > 1, ]
  
  cat(sprintf("  %-25s | Total: %5d | Significant: %d\n",
              comp, nrow(res), nrow(sig)))
  
  write.csv(res, file.path(output_dir, "LIMMA_Results", paste0("DEG_", comp, ".csv")), row.names = FALSE)
  write.csv(sig, file.path(output_dir, "LIMMA_Results", paste0("DEG_SIGNIFICANT_", comp, ".csv")), row.names = FALSE)
}





saveRDS(expr_train, file.path(output_dir, "expr_train.rds"))
saveRDS(pheno_train, file.path(output_dir, "pheno_train.rds"))

saveRDS(expr_test, file.path(output_dir, "expr_test.rds"))
saveRDS(pheno_test, file.path(output_dir, "pheno_test.rds"))

saveRDS(expr_clean, file.path(output_dir, "expr_clean.rds"))      # log2 + filtered
saveRDS(pheno_clean, file.path(output_dir, "pheno_clean.rds"))    # filtered pheno



#Figure 

library(ggplot2)
volcano_plots <- list()
titles <- c("Adenoma_vs_Normal" = "A. Adenoma vs Normal",
            "Cancer_vs_Adenoma" = "B. Cancer vs Adenoma",
            "Metastasis_vs_Cancer" = "C. Metastasis vs Cancer")

for (comp in colnames(contrast.matrix)) {
  res <- topTable(fit2, coef = comp, number = Inf, adjust.method = "BH")
  res$Sig <- ifelse(res$adj.P.Val < 0.05 & abs(res$logFC) > 1,
                    ifelse(res$logFC > 0, "Up", "Down"), "NS")
  p <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = Sig)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c(Up = "red", Down = "blue", NS = "grey70")) +
    theme_minimal(base_size = 13) +
    labs(title = titles[comp], x = "log2 Fold Change", y = "-log10(adj. P-value)") +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))
  volcano_plots[[comp]] <- p
  ggsave(file.path(output_dir, paste0("Figure2_", comp, ".png")), p, width = 6, height = 5, dpi = 300)
}

#Figure End



