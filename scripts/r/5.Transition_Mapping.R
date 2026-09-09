
# STEP INDEX TO STAGE CONVERTER 

library(ggplot2) 

ordered_samples <- readRDS("ordered_samples.rds")

pheno_clean <- readRDS("pheno_clean.rds")





step_results <- read.csv("StepMiner_Q1_Final_Results.csv")

# Stage levels 
stage_order <- c(
  "Normal Colon",
  "Adenoma (Polyp)",
  "Primary Colorectal Cancer",
  "Metastasis (Liver/Lung)"
)

pheno_clean$Stage <- factor(pheno_clean$Stage, levels = stage_order)






# sample-to-stage mapping

sample_stage_map <- data.frame(
  SampleOrder = 1:length(ordered_samples),
  SampleID = ordered_samples,
  Stage = pheno_clean[ordered_samples, "Stage"]
)


# Find stage boundaries 

stage_boundaries <- data.frame()
stage_levels <- levels(pheno_clean$Stage)

for (i in 1:(length(stage_levels) - 1)) {
  
  current_stage <- stage_levels[i]
  next_stage <- stage_levels[i + 1]
  
  # Find last sample of current stage
  last_idx <- max(sample_stage_map$SampleOrder[sample_stage_map$Stage == current_stage])
  
  stage_boundaries <- rbind(stage_boundaries, data.frame(
    Boundary_Index = last_idx,
    From_Stage = current_stage,
    To_Stage = next_stage,
    Transition = paste(current_stage, "→", next_stage)
  ))
}

cat("\n Stage boundaries detected:\n")
print(stage_boundaries)


# Convert each StepIndex to Transition

assign_transition <- function(step_idx, boundaries) {
  
  if (is.na(step_idx)) return(NA)
  
  # Find which boundary this step belongs to
  for (i in 1:nrow(boundaries)) {
    if (step_idx <= boundaries$Boundary_Index[i]) {
      return(as.character(boundaries$Transition[i]))
    }
  }
  
  # If step_idx is after last boundary 
  return(as.character(boundaries$Transition[nrow(boundaries)]))
}


# Apply conversion

step_results$Transition <- sapply(step_results$StepIndex, assign_transition, boundaries = stage_boundaries)


# Also add which stage the step occurs AT (the "before" stage)

get_stage_at_step <- function(step_idx, sample_map) {
  if (is.na(step_idx)) return(NA)
  stage_at_idx <- sample_map$Stage[sample_map$SampleOrder == step_idx]
  if (length(stage_at_idx) == 0) return(NA)
  return(as.character(stage_at_idx))
}

step_results$Step_Occurs_At_Stage <- sapply(step_results$StepIndex, get_stage_at_step, sample_map = sample_stage_map)


# Summary statistics by transition


transition_counts <- table(step_results$Transition)
print(transition_counts)
print(table(step_results$Transition, step_results$Direction))


# Reorder columns for clarity

step_results <- step_results[, c("ProbeID", "SYMBOL", "Transition", "Step_Occurs_At_Stage", 
                                 "StepIndex", "Direction", "F_Stat")]


# Save updated file with step index information

write.csv(step_results, "StepMiner_Q8_Final_Results.csv", row.names = FALSE)





saveRDS(stage_boundaries, "stage_boundaries.rds")
saveRDS(ordered_samples, "ordered_samples.rds")


saveRDS(pheno_clean, "pheno_clean.rds")



#Figure 
trans_df <- as.data.frame(table(step_results$Transition, step_results$Direction))
colnames(trans_df) <- c("Transition", "Direction", "Count")
p3 <- ggplot(trans_df, aes(x = Transition, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal(base_size = 12) +
  coord_flip() +
  labs(title = "CRC Progression Transition (Leakage-Free)", x = "", y = "Gene Count") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
print(p3)
ggsave("Figure3_CRC_Progression_Transition.png", p3, width = 9, height = 5, dpi = 300)

#Figure End

