#LIVER/LUNG MARKER OVERLAP CHECK
sink("Liver_Lung_Overlap_Check.txt")

step_results <- read.csv("StepMiner_Q8_Final_Results.csv")
liver_markers <- c("ALB", "APOA1", "APOB", "TTR", "CYP3A4", "HNF4A", "SERPINA1")
lung_markers  <- c("SFTPC", "SFTPB", "SFTPA1", "NAPSA", "SCGB1A1")
metastasis_genes <- step_results$SYMBOL[
  step_results$Transition == "Primary Colorectal Cancer → Metastasis (Liver/Lung)"
]
metastasis_genes <- unique(metastasis_genes[!is.na(metastasis_genes)])
liver_overlap <- intersect(metastasis_genes, liver_markers)
lung_overlap  <- intersect(metastasis_genes, lung_markers)
cat("Cancer→Metastasis transition genes:", length(metastasis_genes), "\n")
cat("Overlap with liver markers:", length(liver_overlap), "-", paste(liver_overlap, collapse=", "), "\n")
cat("Overlap with lung markers:", length(lung_overlap), "-", paste(lung_overlap, collapse=", "), "\n")

sink()
