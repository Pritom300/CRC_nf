


# PAIRED PRIMARY-VS-METASTASIS CHECK



rm(list = ls())

library(GEOquery)






gse_file <- "GSE41258_series_matrix.txt.gz"

if (!file.exists(gse_file)) {
  stop("File not found: ", gse_file)
}



gse <- getGEO(
  filename = gse_file,
  getGPL = FALSE
)





# LOAD STEPMINER RESULTS


step_results <- read.csv(
  "StepMiner_Q8_Final_Results.csv",
  stringsAsFactors = FALSE
)




raw_pheno <- pData(gse)

#Detecting tissue column.

tissue_col <- names(raw_pheno)[
  sapply(
    raw_pheno,
    function(x) any(grepl("^tissue:", x))
  )
][1]

if (is.na(tissue_col) || length(tissue_col) == 0) {
  stop("Tissue column could not be detected!")
}

cat("Tissue column:", tissue_col, "\n")

raw_pheno$Tissue_Detail <- sub(
  "^tissue: ",
  "",
  raw_pheno[[tissue_col]]
)

cat("\nTissue types:\n")
print(table(raw_pheno$Tissue_Detail))



# PATIENT IDs


patient_col <- names(raw_pheno)[
  sapply(
    raw_pheno,
    function(x) any(grepl("^patient id:", x))
  )
][1]

if (is.na(patient_col) || length(patient_col) == 0) {
  stop("Patient ID column could not be detected!")
}

cat("\nPatient column:", patient_col, "\n")

raw_pheno$Patient_ID <- sub(
  "^patient id: ",
  "",
  raw_pheno[[patient_col]]
)



# METASTASIS PROBES


metastasis_probes <- step_results$ProbeID[
  step_results$Transition ==
    "Primary Colorectal Cancer → Metastasis (Liver/Lung)"
]

metastasis_probes <- unique(
  metastasis_probes[!is.na(metastasis_probes)]
)

cat(
  "\nMetastasis probes from StepMiner:",
  length(metastasis_probes),
  "\n"
)



# FULL EXPRESSION MATRIX


expr_full <- exprs(gse)

if (max(expr_full, na.rm = TRUE) > 50) {
  
  cat("Expression appears unlogged — applying log2(x + 1)\n")
  
  expr_full <- log2(expr_full + 1)
  
} else {
  
  cat("Expression already appears log2-scaled\n")
  
}

cat(
  "Expression matrix:",
  nrow(expr_full),
  "probes x",
  ncol(expr_full),
  "samples\n"
)



# MATCHED PATIENTS


matched_patients <- c(
  "C0136",
  "C0353",
  "C0415",
  "C0575",
  "C0667"
)



# PAIRED PRIMARY vs LIVER METASTASIS (don't have any paired primary and LUNG)


available_probes <- intersect(
  metastasis_probes,
  rownames(expr_full)
)

cat(
  "\nAvailable metastasis probes:",
  length(available_probes),
  "\n"
)

if (length(available_probes) == 0) {
  stop("No metastasis probes found!")
}


pair_results_list <- list()


for (pid in matched_patients) {
  
  cat("\nProcessing patient:", pid, "\n")
  
  primary_id <- rownames(raw_pheno)[
    raw_pheno$Patient_ID == pid &
      raw_pheno$Tissue_Detail == "Primary Tumor"
  ]
  
  met_id <- rownames(raw_pheno)[
    raw_pheno$Patient_ID == pid &
      raw_pheno$Tissue_Detail == "Liver Metastasis"
  ]
  
  
  if (length(primary_id) == 0) {
    
    cat("Primary sample not found\n")
    next
    
  }
  
  
  if (length(met_id) == 0) {
    
    cat("Liver metastasis sample not found\n")
    next
    
  }
  
  
  primary_id <- primary_id[1]
  met_id <- met_id[1]
  
  
  cat("   Primary:", primary_id, "\n")
  cat("   Metastasis:", met_id, "\n")
  
  
  # Metastasis - Primary
  diff_vec <-
    expr_full[available_probes, met_id] -
    expr_full[available_probes, primary_id]
  
  
  pair_results_list[[pid]] <- data.frame(
    
    Patient = pid,
    
    ProbeID = names(diff_vec),
    
    Log2FC_Met_vs_Primary = as.numeric(diff_vec),
    
    stringsAsFactors = FALSE
    
  )
}



# COMBINE RESULTS


if (length(pair_results_list) == 0) {
  
  stop(
    "No patient pairs were successfully found.\n",
    "Check Patient_ID and Tissue_Detail labels."
  )
  
}


pair_results <- do.call(
  rbind,
  pair_results_list
)


cat(
  "\nPaired comparisons:",
  nrow(pair_results),
  "\n"
)


cat("\npair_results columns:\n")

print(colnames(pair_results))


# Critical safety check
if (!"Log2FC_Met_vs_Primary" %in% colnames(pair_results)) {
  
  stop("Log2FC_Met_vs_Primary column was not created!")
  
}



#  SYMBOL + DISCOVERY DIRECTION


pair_results$SYMBOL <-
  
  step_results$SYMBOL[
    match(
      pair_results$ProbeID,
      step_results$ProbeID
    )
  ]


pair_results$Discovery_Direction <-
  
  step_results$Direction[
    match(
      pair_results$ProbeID,
      step_results$ProbeID
    )
  ]






pair_results$Agrees_With_Discovery <-
  
  (
    pair_results$Discovery_Direction == "Up-Step" &
      
      pair_results$Log2FC_Met_vs_Primary > 0
  ) |
  
  (
    pair_results$Discovery_Direction == "Down-Step" &
      
      pair_results$Log2FC_Met_vs_Primary < 0
  )




write.csv(
  pair_results,
  "Metastasis_Genes_Paired_Patient_Check.csv",
  row.names = FALSE
)



# GENE-WISE AGREEMENT


agreement_summary <- aggregate(
  
  Agrees_With_Discovery ~ SYMBOL,
  
  data = pair_results,
  
  FUN = mean
  
)


agreement_summary$Agreement_Percent <-
  
  agreement_summary$Agrees_With_Discovery * 100


agreement_summary <-
  
  agreement_summary[
    order(-agreement_summary$Agreement_Percent),
  ]


write.csv(
  agreement_summary,
  "Metastasis_Genes_Paired_Agreement_Summary.csv",
  row.names = FALSE
)



# FINAL 


sink("Paired_Patient_Check_Summary.txt")

overall_agreement <-
  
  mean(
    pair_results$Agrees_With_Discovery,
    na.rm = TRUE
  ) * 100


cat("\n", paste(rep("=", 70), collapse = ""), "\n")

cat(
  "\nMean agreement rate across",
  length(unique(pair_results$Patient)),
  "matched patients:\n"
)

cat(
  "(fraction of patients where within-patient direction matches discovery direction)\n"
)

cat(
  "Overall mean agreement:",
  round(overall_agreement, 1),
  "%\n"
)

cat(
  "\nPAIRED PRIMARY-vs-METASTASIS CHECK COMPLETE!\n"
)




#Additional Check

table(raw_pheno$Tissue_Detail)





# only Primary / Liver / Lung samples
pair_check <- raw_pheno[
  raw_pheno$Tissue_Detail %in% c(
    "Primary Tumor",
    "Liver Metastasis",
    "Lung Metastasis"
  ),
  c("Patient_ID", "Tissue_Detail")
]

# per-patient tissue count
patient_tissues <- aggregate(
  Tissue_Detail ~ Patient_ID,
  data = pair_check,
  FUN = function(x) paste(unique(x), collapse = " | ")
)

print(patient_tissues)


# how many patients have Primary + Liver?
primary_liver_patients <- unique(
  pair_check$Patient_ID[
    pair_check$Patient_ID %in%
      pair_check$Patient_ID[pair_check$Tissue_Detail == "Primary Tumor"] &
      pair_check$Patient_ID %in%
      pair_check$Patient_ID[pair_check$Tissue_Detail == "Liver Metastasis"]
  ]
)

# Checking how many patients have Primary + Lung
primary_lung_patients <- unique(
  pair_check$Patient_ID[
    pair_check$Patient_ID %in%
      pair_check$Patient_ID[pair_check$Tissue_Detail == "Primary Tumor"] &
      pair_check$Patient_ID %in%
      pair_check$Patient_ID[pair_check$Tissue_Detail == "Lung Metastasis"]
  ]
)

cat("Primary + Liver matched patients:", 
    length(primary_liver_patients), "\n")

cat("Primary + Lung matched patients:", 
    length(primary_lung_patients), "\n")

cat("\nPrimary + Liver IDs:\n")
print(primary_liver_patients)

cat("\nPrimary + Lung IDs:\n")
print(primary_lung_patients)

sink()
