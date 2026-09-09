# Probe Annotation and Mapping

library(AnnotationDbi)
library(hgu133plus2.db)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE) 
input_file <- args[1]
output_dir <- args[2]

map_probes_q1 <- function(file_path, out_dir) {
  

  
  if (!file.exists(file_path)) {
    cat("File not found\n")
    return(NULL)
  }
  
  deg <- read.csv(file_path, stringsAsFactors = FALSE)
  
  if(!"ProbeID" %in% colnames(deg)){
    stop("ProbeID column missing!")
  }
  
  deg$ProbeID <- as.character(deg$ProbeID)
  
  # Mapping Started
  mapping <- AnnotationDbi::select(
    hgu133plus2.db,
    keys = deg$ProbeID,
    columns = c("SYMBOL", "GENENAME"),
    keytype = "PROBEID"
  )
  
  # Merge
  deg_mapped <- merge(deg, mapping,
                      by.x = "ProbeID",
                      by.y = "PROBEID",
                      all.x = TRUE)
  
  # REMOVING NA Genes
  deg_mapped <- deg_mapped[!is.na(deg_mapped$SYMBOL), ]
  
  # HANDLE MULTIPLE PROBES PER GENE
  deg_mapped <- deg_mapped %>%
    group_by(SYMBOL) %>%
    slice_max(order_by = abs(logFC), n = 1) %>%
    ungroup() 
  
  # Sort by significance
  deg_mapped <- deg_mapped[order(deg_mapped$adj.P.Val), ]
  
  base_name <- basename(file_path)
  new_name <- sub("\\.csv$", "_Q1Mapped.csv", base_name)
  
  write.csv(deg_mapped, file.path(out_dir, new_name), row.names = FALSE)
  
  cat("Final genes:", nrow(deg_mapped), "\n")
  cat("Unique genes:", length(unique(deg_mapped$SYMBOL)), "\n")
  
  return(deg_mapped)
}



if (!grepl("SIGNIFICANT", input_file)) {
  map_probes_q1(input_file, output_dir)
}

cat("\n DONE! gene list created!\n")
