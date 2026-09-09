# BOOTSTRAP 95% CONFIDENCE INTERVALS

library(randomForest)
set.seed(123)

n_boot_ci <- 2000

rf_model <- readRDS("Stage_Classification_Final.rds")
expr_test <- readRDS("expr_test.rds")
pheno_test <- readRDS("pheno_test.rds")
common_genes <- readRDS("common_genes.rds")
X_test <- readRDS("X_test.rds")
y_test <- readRDS("y_test.rds")
predictions <- readRDS("predictions_GSE41333.rds")
accuracy <- readRDS("accuracy.rds")
accuracy_primary <- readRDS("accuracy_primary.rds")


bootstrap_ci <- function(true_labels, pred_labels, n_boot = n_boot_ci) {
  n <- length(true_labels)
  acc_boot <- numeric(n_boot)
  
  classes <- levels(as.factor(true_labels))
  sens_boot <- matrix(NA, n_boot, length(classes), dimnames = list(NULL, classes))
  prec_boot <- matrix(NA, n_boot, length(classes), dimnames = list(NULL, classes))
  f1_boot   <- matrix(NA, n_boot, length(classes), dimnames = list(NULL, classes))
  
  for (i in seq_len(n_boot)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    t_b <- true_labels[idx]; p_b <- pred_labels[idx]
    acc_boot[i] <- mean(t_b == p_b)
    
    for (c in classes) {
      tp <- sum(p_b == c & t_b == c)
      fp <- sum(p_b == c & t_b != c)
      fn <- sum(p_b != c & t_b == c)
      sens <- ifelse(tp + fn > 0, tp / (tp + fn), NA)
      prec <- ifelse(tp + fp > 0, tp / (tp + fp), NA)
      f1   <- ifelse(!is.na(sens) && !is.na(prec) && (sens + prec) > 0,
                     2 * sens * prec / (sens + prec), NA)
      sens_boot[i, c] <- sens; prec_boot[i, c] <- prec; f1_boot[i, c] <- f1
    }
  }
  
  ci <- function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE)
  
  list(
    Accuracy_CI = ci(acc_boot),
    PerClass = data.frame(
      Class = classes,
      Sensitivity_Lower = sapply(classes, function(c) ci(sens_boot[, c])[1]),
      Sensitivity_Upper = sapply(classes, function(c) ci(sens_boot[, c])[2]),
      Precision_Lower   = sapply(classes, function(c) ci(prec_boot[, c])[1]),
      Precision_Upper   = sapply(classes, function(c) ci(prec_boot[, c])[2]),
      F1_Lower          = sapply(classes, function(c) ci(f1_boot[, c])[1]),
      F1_Upper          = sapply(classes, function(c) ci(f1_boot[, c])[2])
    )
  )
}


# HELD-OUT TEST SET — 95% Bootstrap CI

predictions_test <- predict(rf_model, newdata = X_test)
ci_test <- bootstrap_ci(as.character(y_test), as.character(predictions_test))

cat("Overall accuracy:", round(accuracy, 3),
    " 95% CI: [", round(ci_test$Accuracy_CI[1], 3), ",", round(ci_test$Accuracy_CI[2], 3), "]\n")
    
print(ci_test$PerClass)

write.csv(ci_test$PerClass, "Test_Set_Bootstrap_CI.csv", row.names = FALSE)

#  GSE14333 external validation — 95% Bootstrap CI 

n_ext <- length(predictions)  
acc_ext_boot <- numeric(n_boot_ci)

for (i in seq_len(n_boot_ci)) {
  idx <- sample(seq_len(n_ext), size = n_ext, replace = TRUE)
  acc_ext_boot[i] <- mean(predictions[idx] == "Primary Colorectal Cancer")
}

ci_ext <- quantile(acc_ext_boot, c(0.025, 0.975))
cat("GSE14333 accuracy:", round(accuracy_primary, 3),
    " 95% CI: [", round(ci_ext[1], 3), ",", round(ci_ext[2], 3), "]\n")

write.csv(
  data.frame(
    Accuracy = accuracy_primary,
    CI_Lower = ci_ext[1],
    CI_Upper = ci_ext[2]
  ),
  "External_Validation_Bootstrap_CI.csv",
  row.names = FALSE
)
