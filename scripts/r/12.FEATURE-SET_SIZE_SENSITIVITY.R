# FEATURE-SET SIZE SENSITIVITY (tuned within TRAIN only, via 5-fold CV)
library(caret)
library(ggplot2)
library(randomForest)
set.seed(123)


pheno_clean <- readRDS("pheno_clean.rds")
expr_clean <- readRDS("expr_clean.rds")
step_results <- read.csv(
  "StepMiner_Q8_Final_Results.csv")
  
candidate_n <- c(10, 15, 20, 25, 30, 40, 50)

# 5-fold CV folds within TRAIN samples only (test set never touched here)

cv_folds <- createFolds(pheno_clean$Stage, k = 5, list = TRUE)


get_top_genes_n <- function(transition_name, step_df, n) {
  genes <- step_df$SYMBOL[step_df$Transition == transition_name]
  genes <- unique(genes[!is.na(genes)])
  gene_fstat <- step_df[step_df$SYMBOL %in% genes, c("SYMBOL", "F_Stat")]
  gene_fstat <- gene_fstat[!duplicated(gene_fstat$SYMBOL), ]
  gene_fstat <- gene_fstat[order(-gene_fstat$F_Stat), ]
  head(gene_fstat$SYMBOL, n)
}
sensitivity_results <- data.frame()
for (n in candidate_n) {
  
  selected_n <- unique(c(
    get_top_genes_n("Normal Colon → Adenoma (Polyp)", step_results, n),
    get_top_genes_n("Adenoma (Polyp) → Primary Colorectal Cancer", step_results, n),
    get_top_genes_n("Primary Colorectal Cancer → Metastasis (Liver/Lung)", step_results, n)
  ))
  selected_n <- intersect(selected_n, rownames(expr_clean))
  
  fold_acc <- c()
  
  for (fold_name in names(cv_folds)) {
    val_idx   <- cv_folds[[fold_name]]
    val_samples   <- rownames(pheno_clean)[val_idx]
    train_samples_cv <- setdiff(rownames(pheno_clean), val_samples)
    
    Xc_train <- as.data.frame(t(expr_clean[selected_n, train_samples_cv, drop = FALSE]))
    yc_train <- droplevels(as.factor(pheno_clean[train_samples_cv, "Stage"]))
    Xc_val   <- as.data.frame(t(expr_clean[selected_n, val_samples, drop = FALSE]))
    yc_val   <- droplevels(as.factor(pheno_clean[val_samples, "Stage"]))
    colnames(Xc_train) <- make.names(colnames(Xc_train))
    colnames(Xc_val)   <- make.names(colnames(Xc_val))
    
    m <- randomForest(x = Xc_train, y = yc_train, ntree = 500)
    pred <- predict(m, newdata = Xc_val)
    fold_acc <- c(fold_acc, mean(pred == yc_val))
  }
  
  sensitivity_results <- rbind(sensitivity_results, data.frame(
    n_per_transition = n,
    n_unique_genes = length(selected_n),
    Mean_CV_Accuracy = mean(fold_acc),
    SD_CV_Accuracy = sd(fold_acc)
  ))
  
  cat("n =", n, "-> unique genes:", length(selected_n),
      "| mean CV accuracy:", round(mean(fold_acc), 3),
      "± ", round(sd(fold_acc), 3), "\n")
}


write.csv(sensitivity_results, "Feature_Size_Sensitivity.csv", row.names = FALSE)
cat("\nFEATURE-SIZE SENSITIVITY SUMMARY:\n")
print(sensitivity_results)
best_row <- sensitivity_results[which.max(sensitivity_results$Mean_CV_Accuracy), ]
cat("\nBest n_per_transition by CV accuracy:", best_row$n_per_transition,
    "(", best_row$n_unique_genes, "unique genes )\n")
    
    
    
#Figure 

p_s2 <- ggplot(sensitivity_results, aes(x = n_per_transition, y = Mean_CV_Accuracy)) +
  geom_line(color = "darkblue") +
  geom_point(size = 2, color = "darkblue") +
  geom_errorbar(aes(ymin = Mean_CV_Accuracy - SD_CV_Accuracy, ymax = Mean_CV_Accuracy + SD_CV_Accuracy), width = 1) +
  geom_vline(xintercept = 25, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 12) +
  labs(title = "Feature-Set Size Sensitivity (5-fold CV, train only)",
       x = "Genes Selected per Transition", y = "Mean CV Accuracy") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
print(p_s2)
ggsave("FigureS2_Feature_Size_Sensitivity.png", p_s2, width = 8, height = 5, dpi = 300)

#Figure End
