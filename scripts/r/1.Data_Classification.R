
if (!requireNamespace("GEOquery", quietly=TRUE)) {
  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
  BiocManager::install("GEOquery")
}
library(GEOquery)



args <- commandArgs(trailingOnly=TRUE)
gse_file <- args[1]
output_dir <- args[2]
if (!file.exists(gse_file)) {
  stop("File not found!")
}




# Load data 

gse <- tryCatch(
  getGEO(filename=gse_file, getGPL=FALSE),
  error = function(e) stop("Failed to load! ", e$message)
)

pheno <- pData(gse)


# Extract metadata 

meta_cols <- grep("title|characteristics", colnames(pheno), value=TRUE, ignore.case=TRUE)


cat(paste(meta_cols, collapse=", "), "\n")

combined_text <- apply(pheno[, meta_cols, drop=FALSE], 1, function(x) {
  x[is.na(x)] <- ""   
  paste(tolower(trimws(x)), collapse=" ")
})


pattern_meta <- "\\b(metastasis|metastatic|liver met|hepatic metastasis|lung met|pulmonary metastasis|secondary)\\b"
pattern_primary <- "\\b(primary tumor|primary colorectal cancer|crc|colorectal carcinoma|adenocarcinoma)\\b"
pattern_adenoma <- "\\b(adenoma|polyp|tubular|villous|serrated|hyperplastic)\\b"
pattern_normal <- "\\b(normal|healthy|control|non[- ]neoplastic|adjacent mucosa|uninvolved)\\b"


# Classification (priority-based)

pheno$Stage <- "Unclassified/Other"

pheno$Stage[grepl(pattern_meta, combined_text)] <- "Metastasis (Liver/Lung)"

pheno$Stage[grepl(pattern_primary, combined_text) & pheno$Stage=="Unclassified/Other"] <- "Primary Colorectal Cancer"

pheno$Stage[grepl(pattern_adenoma, combined_text) & pheno$Stage=="Unclassified/Other"] <- "Adenoma (Polyp)"

pheno$Stage[grepl(pattern_normal, combined_text) & pheno$Stage=="Unclassified/Other"] <- "Normal Colon"


# ordering

stage_order <- c(
  "Normal Colon",
  "Adenoma (Polyp)",
  "Primary Colorectal Cancer",
  "Metastasis (Liver/Lung)",
  "Unclassified/Other"
)

pheno$Stage <- factor(pheno$Stage, levels=stage_order)


unclassified_ids <- rownames(pheno)[pheno$Stage == "Unclassified/Other"]



if(length(unclassified_ids) > 0){
  write.csv(
   data.frame(SampleID = unclassified_ids),
    file.path(output_dir, "Unclassified_Samples.csv"),
    row.names = FALSE
  )
}


counts <- table(pheno$Stage)
total_samples <- sum(counts)
percentages <- round(100 * counts / total_samples, 1)

summary_df <- data.frame(
  Stage = names(counts),
  Count = as.numeric(counts),
  Percent = percentages
)




pub_table <- data.frame(
  Stage = stage_order,
  Description = c(
    "Healthy colonic mucosa",
    "Pre-malignant adenoma",
    "Primary colorectal cancer",
    "Metastatic lesions",
    "Ambiguous / excluded samples"
  ),
  n = as.numeric(counts),
  Percent = percentages,
  check.names = FALSE
)

cat("\n Publication Table:\n")
print(pub_table, row.names = FALSE)


write.csv(pub_table, file.path(output_dir, "GSE41258_Sample_Table.csv"), row.names = FALSE)
write.csv(pheno, file.path(output_dir, "GSE41258_Classified_Metadata.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))







