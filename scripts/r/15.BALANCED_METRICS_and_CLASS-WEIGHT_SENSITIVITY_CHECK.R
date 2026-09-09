#  BALANCED METRICS + CLASS-WEIGHT SENSITIVITY CHECK

library(pROC)
library(randomForest)
library(pROC)


X_train <- readRDS("X_train.rds")
y_train <- readRDS("y_train.rds")
X_test <- readRDS("X_test.rds")
y_test <- readRDS("y_test.rds")
rf_model <- readRDS("Stage_Classification_Final.rds")

predictions <- predict(rf_model, newdata = X_test)


confusion <- table(
  Predicted = predictions,
  Actual = y_test
)
accuracy <- sum(diag(confusion)) / sum(confusion)
classes <- levels(y_test)

sink("Balanced_Metrics_Class_Weight_Summary.txt")

#  Balanced accuracy & macro-F1 from the EXISTING (unweighted) model 

per_class_sens <- sapply(classes, function(c) {
  tp <- if(c %in% rownames(confusion) & c %in% colnames(confusion)) confusion[c, c] else 0
  fn <- if(c %in% colnames(confusion)) sum(confusion[, c]) - tp else 0
  ifelse(tp + fn > 0, tp / (tp + fn), 0)
})

per_class_prec <- sapply(classes, function(c) {
  tp <- if(c %in% rownames(confusion) & c %in% colnames(confusion)) confusion[c, c] else 0
  fp <- if(c %in% rownames(confusion)) sum(confusion[c, ]) - tp else 0
  ifelse(tp + fp > 0, tp / (tp + fp), 0)
})

per_class_f1 <- ifelse((per_class_sens + per_class_prec) > 0,
                       2 * per_class_sens * per_class_prec / (per_class_sens + per_class_prec), 0)
                       
balanced_accuracy <- mean(per_class_sens)

macro_f1 <- mean(per_class_f1)

cat("UNWEIGHTED MODEL (current)\n")
cat("Overall accuracy:  ", round(accuracy, 3), "\n")
cat("Balanced accuracy: ", round(balanced_accuracy, 3), "\n")
cat("Macro-F1:          ", round(macro_f1, 3), "\n")



# Multiclass AUC (one-vs-rest) 

test_probs <- predict(rf_model, newdata = X_test, type = "prob")
mc_auc <- multiclass.roc(y_test, test_probs)

cat("Multiclass AUC (one-vs-rest, mean): ", round(as.numeric(mc_auc$auc), 3), "\n\n")

# SENSITIVITY CHECK:  (does class-weighting change things or not)
# Inverse-frequency class weights, computed from TRAIN labels only 

class_counts <- table(y_train)
class_weights <- max(class_counts) / class_counts   # inverse-frequency



print(round(class_weights, 2))

set.seed(123)

rf_weighted <- randomForest(
  x = X_train, y = y_train,
  ntree = 500, importance = TRUE,
  classwt = as.numeric(class_weights)
)


pred_weighted <- predict(rf_weighted, newdata = X_test)
confusion_w <- table(Predicted = pred_weighted, Actual = y_test)
sens_w <- sapply(classes, function(c) {
  tp <- if(c %in% rownames(confusion_w) & c %in% colnames(confusion_w)) confusion_w[c, c] else 0
  fn <- if(c %in% colnames(confusion_w)) sum(confusion_w[, c]) - tp else 0
  ifelse(tp + fn > 0, tp / (tp + fn), 0)
})


prec_w <- sapply(classes, function(c) {
  tp <- if(c %in% rownames(confusion_w) & c %in% colnames(confusion_w)) confusion_w[c, c] else 0
  fp <- if(c %in% rownames(confusion_w)) sum(confusion_w[c, ]) - tp else 0
  ifelse(tp + fp > 0, tp / (tp + fp), 0)
})


f1_w <- ifelse((sens_w + prec_w) > 0, 2 * sens_w * prec_w / (sens_w + prec_w), 0)

accuracy_w <- sum(diag(confusion_w)) / sum(confusion_w)

balanced_accuracy_w <- mean(sens_w)

macro_f1_w <- mean(f1_w)


cat("\nCLASS-WEIGHTED MODEL (sensitivity check) \n")
cat("Overall accuracy:  ", round(accuracy_w, 3), "\n")
cat("Balanced accuracy: ", round(balanced_accuracy_w, 3), "\n")
cat("Macro-F1:          ", round(macro_f1_w, 3), "\n")
cat("Per-class sensitivity (weighted):\n")
print(round(sens_w, 3))

sink()

# side-by-side comparison table 
comparison <- data.frame(
  Class = classes,
  Sensitivity_Unweighted = round(per_class_sens, 3),
  Sensitivity_Weighted   = round(sens_w, 3),
  F1_Unweighted = round(per_class_f1, 3),
  F1_Weighted   = round(f1_w, 3)
)

write.csv(comparison, "Class_Weighting_Comparison.csv", row.names = FALSE)
print(comparison)
