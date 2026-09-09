
#WITHIN-STAGE ORDER SENSITIVITY / PERMUTATION CHECK
library(ggplot2)




pheno_clean <- readRDS("pheno_clean.rds") 
ordered_samples <- readRDS("ordered_samples.rds")
expr_final <- readRDS("expr_final.rds") 
stage_boundaries <- readRDS("stage_boundaries.rds")  
step_results <- readRDS("step_results.rds") 
selected_genes <- readRDS("selected_genes.rds")

set.seed(123)
n_perm <- 100

# vectorized version of fitStep, we used here because we need ~150,000 calls

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

# Shuffle sample order WITHIN each stage block
permute_within_stage <- function(stage_vec) {
  idx <- seq_along(stage_vec)
  for (s in unique(stage_vec)) {
    block <- which(stage_vec == s)
    idx[block] <- sample(block)
  }
  idx
}


# Convert each StepIndex to Transition

assign_transition <- function(step_idx, boundaries) {
  
  if (is.na(step_idx)) return(NA)
  
  # Find which boundary this step belongs to
  for (i in 1:nrow(boundaries)) {
    if (step_idx <= boundaries$Boundary_Index[i]) {
      return(as.character(boundaries$Transition[i]))
    }
  }
  

  return(as.character(boundaries$Transition[nrow(boundaries)]))
}



stage_vec <- as.character(pheno_clean[ordered_samples, "Stage"])
genes_to_test <- rownames(expr_final)   # all monotonic genes



perm_results <- vector("list", length(genes_to_test))

for (k in seq_along(genes_to_test)) {
  gene <- genes_to_test[k]
  vec <- as.numeric(expr_final[gene, ])
  
  obs <- fitStep_fast(vec)
  obs_transition <- assign_transition(obs$index, stage_boundaries)
  
  perm_idx <- integer(n_perm)
  for (p in seq_len(n_perm)) {
    ord <- permute_within_stage(stage_vec)
    perm_idx[p] <- fitStep_fast(vec[ord])$index
  }
  
  perm_transitions <- sapply(perm_idx, assign_transition, boundaries = stage_boundaries)
  
  perm_results[[k]] <- data.frame(
    ProbeID = gene,
    Observed_StepIndex = obs$index,
    Observed_Transition = obs_transition,
    Pct_Same_Transition = mean(perm_transitions == obs_transition),
    StepIndex_SD = sd(perm_idx)
  )
  
  if (k %% 200 == 0) cat("   processed", k, "/", length(genes_to_test), "genes\n")
}

perm_results <- do.call(rbind, perm_results)
write.csv(perm_results, "StepMiner_Permutation_Stability.csv", row.names = FALSE)

cat("\nPERMUTATION STABILITY SUMMARY:\n")
cat("   Genes tested:", nrow(perm_results), "\n")
cat("   Mean % of permutations landing in the SAME transition as observed:",
    round(mean(perm_results$Pct_Same_Transition) * 100, 1), "%\n")
cat("   Genes 100% stable across all", n_perm, "permutations:",
    sum(perm_results$Pct_Same_Transition == 1), "\n")
cat("   Genes where transition changes in >50% of permutations (unstable):",
    sum(perm_results$Pct_Same_Transition < 0.5), "\n")

# Also check specifically for the 72 genes used in the final RF model
final_gene_check <- perm_results[perm_results$ProbeID %in%
                                   step_results$ProbeID[step_results$SYMBOL %in% selected_genes], ]
cat("\n   Among the genes feeding the final biomarker panel — mean stability:",
    round(mean(final_gene_check$Pct_Same_Transition) * 100, 1), "%\n")



# Figure 

p_s1 <- ggplot(perm_results, aes(x = Pct_Same_Transition)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "white") +
  geom_vline(xintercept = mean(perm_results$Pct_Same_Transition), color = "red", linetype = "dashed") +
  theme_minimal(base_size = 12) +
  labs(title = "StepMiner Transition-Assignment Stability (100 within-stage permutations)",
       x = "% of Permutations Matching Observed Transition", y = "Number of Genes") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11))
print(p_s1)
ggsave("FigureS1_StepMiner_Permutation_Stability.png", p_s1, width = 8, height = 5, dpi = 300)

#Figure end
