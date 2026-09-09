# ARE METASTASIS-SIGNATURE GENES ORGAN MARKERS?


if (!requireNamespace("GEOquery", quietly=TRUE)) {
  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
  BiocManager::install("GEOquery")
}
library(GEOquery)

gse_file <- "GSE41258_series_matrix.txt.gz"
if (!file.exists(gse_file)) {
  stop("File not found! ")
}


# Load data 

gse <- tryCatch(
  getGEO(filename=gse_file, getGPL=FALSE),
  error = function(e) stop("Failed to load GEO data: ", e$message)
)
step_results <- read.csv("StepMiner_Q8_Final_Results.csv")
raw_pheno <- pData(gse)
tissue_col <- names(raw_pheno)[
  sapply(raw_pheno, function(x) any(grepl("^tissue:", x)))
][1]
raw_pheno$Tissue_Detail <- sub("^tissue: ", "", raw_pheno[[tissue_col]])
normal_liver_ids <- rownames(raw_pheno)[raw_pheno$Tissue_Detail == "Normal Liver"]
normal_lung_ids  <- rownames(raw_pheno)[raw_pheno$Tissue_Detail == "Normal Lung"]
normal_colon_ids <- rownames(raw_pheno)[raw_pheno$Tissue_Detail == "Normal Colon"]
cat("Normal Liver samples found:", length(normal_liver_ids), "\n")
cat("Normal Lung samples found:", length(normal_lung_ids), "\n")
cat("Normal Colon samples found:", length(normal_colon_ids), "\n")

metastasis_probes <- step_results$ProbeID[
  step_results$Transition == "Primary Colorectal Cancer → Metastasis (Liver/Lung)"
]
expr_full <- exprs(gse)  # raw, all 390 samples, before train/test split
if (max(expr_full, na.rm = TRUE) > 50) expr_full <- log2(expr_full + 1)
get_organ_expr <- function(ids) {
  ids <- intersect(ids, colnames(expr_full))
  rowMeans(expr_full[intersect(metastasis_probes, rownames(expr_full)), ids, drop = FALSE])
}
organ_compare <- data.frame(
  ProbeID = intersect(metastasis_probes, rownames(expr_full)),
  Mean_NormalColon = get_organ_expr(normal_colon_ids),
  Mean_NormalLiver = get_organ_expr(normal_liver_ids),
  Mean_NormalLung  = get_organ_expr(normal_lung_ids)
)
organ_compare$SYMBOL <- step_results$SYMBOL[match(organ_compare$ProbeID, step_results$ProbeID)]
organ_compare$Liver_minus_Colon <- organ_compare$Mean_NormalLiver - organ_compare$Mean_NormalColon
organ_compare$Lung_minus_Colon  <- organ_compare$Mean_NormalLung  - organ_compare$Mean_NormalColon
write.csv(organ_compare, "Metastasis_Genes_Organ_Baseline_Check.csv", row.names = FALSE)

sink("Metastasis_Genes_Organ_Baseline_Summary.txt")
cat("\nGenes with large liver-baseline shift (|diff| > 2 log2):",
    sum(abs(organ_compare$Liver_minus_Colon) > 2, na.rm = TRUE), "/", nrow(organ_compare), "\n")
cat("Genes with large lung-baseline shift (|diff| > 2 log2):",
    sum(abs(organ_compare$Lung_minus_Colon) > 2, na.rm = TRUE), "/", nrow(organ_compare), "\n")
cat("\nTop 10 genes most liver-baseline-shifted (potential organ-identity confound):\n")
print(head(organ_compare[order(-abs(organ_compare$Liver_minus_Colon)), c("SYMBOL","Liver_minus_Colon")], 10))
sink()
