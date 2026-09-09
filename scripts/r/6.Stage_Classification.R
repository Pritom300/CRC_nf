
library(hgu133plus2.db)
library(AnnotationDbi)
library(randomForest)

expr_train <- readRDS("expr_train.rds")
expr_test <- readRDS("expr_test.rds")

pheno_test <- readRDS("pheno_test.rds")

expr_clean <- readRDS("expr_clean.rds")      
pheno_clean <- readRDS("pheno_clean.rds") 

collapse_to_symbol <- function(expr_mat) {
  if (max(expr_mat, na.rm = TRUE) > 50) {
    expr_mat <- log2(expr_mat + 1)
  }
  probe_to_symbol <- mapIds(
    hgu133plus2.db,
    keys = rownames(expr_mat),
    column = "SYMBOL",
    keytype = "PROBEID",
    multiVals = "first"
  )
  keep <- !is.na(probe_to_symbol)
  expr_mat <- expr_mat[keep, ]
  rownames(expr_mat) <- probe_to_symbol[keep]
  expr_mat <- aggregate(expr_mat, by = list(Gene = rownames(expr_mat)), FUN = max)
  rownames(expr_mat) <- expr_mat$Gene
  expr_mat$Gene <- NULL
  as.matrix(expr_mat)
}

# Building symbol-level matrix for TRAIN...

expr_clean <- collapse_to_symbol(expr_train)


# Building symbol-level matrix for TEST

expr_test_sym <- collapse_to_symbol(expr_test)



# STAGE CLASSIFICATION, TRAIN ON TRAIN, EVALUATE ONCE ON HELD-OUT TEST



step_results <- read.csv("StepMiner_Q8_Final_Results.csv")

get_top_genes <- function(transition_name, step_df, n = 25) {
  genes <- step_df$SYMBOL[step_df$Transition == transition_name]
  genes <- unique(genes[!is.na(genes)])
  gene_fstat <- step_df[step_df$SYMBOL %in% genes, c("SYMBOL", "F_Stat")]
  gene_fstat <- gene_fstat[!duplicated(gene_fstat$SYMBOL), ]
  gene_fstat <- gene_fstat[order(-gene_fstat$F_Stat), ]
  return(head(gene_fstat$SYMBOL, n))
}

selected_genes <- unique(c(
  get_top_genes("Normal Colon → Adenoma (Polyp)", step_results, 25),
  get_top_genes("Adenoma (Polyp) → Primary Colorectal Cancer", step_results, 25),
  get_top_genes("Primary Colorectal Cancer → Metastasis (Liver/Lung)", step_results, 25)
))
selected_genes <- selected_genes[!is.na(selected_genes)]  #balanced genes (discovered on TRAIN only)


common_genes <- intersect(selected_genes, intersect(rownames(expr_clean), rownames(expr_test_sym)))  #Common genes available in both train & test:



# TRAIN matrix
X_train <- as.data.frame(t(expr_clean[common_genes, rownames(pheno_clean), drop = FALSE]))
y_train <- droplevels(as.factor(pheno_clean$Stage))
colnames(X_train) <- make.names(colnames(X_train))

# TEST matrix 
test_common_samples <- intersect(colnames(expr_test_sym), rownames(pheno_test))
X_test <- as.data.frame(t(expr_test_sym[common_genes, test_common_samples, drop = FALSE]))
y_test  <- factor(pheno_test[test_common_samples, "Stage"], levels = levels(y_train))
colnames(X_test) <- make.names(colnames(X_test))

cat("\n Train:", nrow(X_train), "| Test:", nrow(X_test), "\n")
cat(" Sample overlap check",
    length(intersect(rownames(pheno_clean), rownames(pheno_test))), "\n")


# TRAIN RANDOM FOREST


set.seed(123)
rf_model <- randomForest(
  x = X_train, y = y_train,
  ntree = 500, importance = TRUE, na.action = na.omit
)




# EVALUATE ONCE ON HELD-OUT TEST


predictions <- predict(rf_model, newdata = X_test)
confusion <- table(Predicted = predictions, Actual = y_test)


print(confusion)

accuracy <- sum(diag(confusion)) / sum(confusion)
cat("\nOVERALL ACCURACY (truly held-out):", round(accuracy * 100, 2), "%\n")

classes <- levels(y_test)
cat("\n PER-CLASS PERFORMANCE:\n")
for(c in classes) {
  tp <- if(c %in% rownames(confusion) & c %in% colnames(confusion)) confusion[c, c] else 0
  fp <- if(c %in% rownames(confusion)) sum(confusion[c, ]) - tp else 0
  fn <- if(c %in% colnames(confusion)) sum(confusion[, c]) - tp else 0
  sens <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
  prec <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
  f1 <- ifelse(sens + prec > 0, 2 * sens * prec / (sens + prec), 0)
  cat(sprintf("   %-25s | Sensitivity: %.3f | Precision: %.3f | F1: %.3f\n", c, sens, prec, f1))
}


# FEATURE IMPORTANCE 


importance_df <- data.frame(
  Gene = rownames(rf_model$importance),
  Gini_Importance = rf_model$importance[, "MeanDecreaseGini"],
  Permutation_Importance = rf_model$importance[, "MeanDecreaseAccuracy"]
)
importance_df <- importance_df[order(-importance_df$Gini_Importance), ]

# Check rank concordance between the two importance measures
importance_df$Gini_Rank <- rank(-importance_df$Gini_Importance)
importance_df$Permutation_Rank <- rank(-importance_df$Permutation_Importance)
rank_correlation <- cor(importance_df$Gini_Rank, importance_df$Permutation_Rank, method = "spearman")

cat("\nTOP 10 IMPORTANT GENES (Gini):\n")
print(head(importance_df[order(-importance_df$Gini_Importance), c("Gene","Gini_Importance","Permutation_Importance")], 10))

cat("\nTOP 10 IMPORTANT GENES (Permutation):\n")
print(head(importance_df[order(-importance_df$Permutation_Importance), c("Gene","Gini_Importance","Permutation_Importance")], 10))

cat("\n Spearman rank correlation between Gini and Permutation importance:",
    round(rank_correlation, 3), "\n")

write.csv(importance_df, "Stage_Classification_Importance.csv", row.names = FALSE)


saveRDS(rf_model, "Stage_Classification_Final.rds")


saveRDS(selected_genes, "selected_genes.rds")
saveRDS(X_train, "X_train.rds")
saveRDS(y_train, "y_train.rds")
saveRDS(X_test, "X_test.rds")
saveRDS(y_test, "y_test.rds")
saveRDS(common_genes, "common_genes.rds") 

saveRDS(expr_train, "expr_train.rds")
saveRDS(expr_test, "expr_test.rds")
saveRDS(pheno_test, "pheno_test.rds")
saveRDS(expr_clean, "expr_clean.rds")
saveRDS(pheno_clean, "pheno_clean.rds")

saveRDS(predictions, "predictions.rds")

saveRDS(accuracy, "accuracy.rds")





#Figure 

# Confusion Matrix Heatmap 


library(ggplot2)
library(reshape2)


conf_mat <- as.data.frame(confusion)
colnames(conf_mat) <- c("Predicted", "Actual", "Freq")


# Plot heatmap

p1 <- ggplot(conf_mat, aes(x = Actual, y = Predicted, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 5, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Stage Classification Confusion Matrix",
    x = "Actual Stage",
    y = "Predicted Stage",
    fill = "Count"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p1)

# Save
ggsave(
  "Figure4A_Confusion_Matrix_Heatmap.png",
  plot = p1,
  width = 10,
  height = 8,
  dpi = 300
)




# Top 15 Feature Importance Plot


importance_df <- read.csv("Stage_Classification_Importance.csv")

top15 <- importance_df[order(-importance_df$Gini_Importance), ][1:15, ]

# Remove duplicate genes 
top15 <- top15[!duplicated(top15$Gene), ]

top15$Gene <- factor(top15$Gene, levels = rev(top15$Gene))


# Plot feature importance

p2 <- ggplot(top15, aes(x = Gene, y = Gini_Importance)) +
  geom_bar(stat = "identity", fill = "darkblue") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Top 15 Stage Classification Biomarkers",
    x = "Gene",
    y = "Mean Decrease Gini"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p2)

# Save
ggsave(
  "Figure4B_Top15_Feature_Importance.png",
  plot = p2,
  width = 10,
  height = 8,
  dpi = 300
)

#Figure end
