library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)
library(caret)
library(caTools)

args <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1])

# Get the parameters from stored Parameters file
parameters <- read.csv("/hpc/shared/prekovic/dhaynessimmons/data/CIN/XGB_cat_parameters.csv")

# select the parameters and weights corrsponding to the index
selected_parameters <- parameters[index, ]
rm(parameters)
cat("\n Parameters: \n")
print(selected_parameters)
cat("\n\n")

selected_feature <- selected_parameters$Feature
selected_rna_set <- selected_parameters$RNA_set
selected_trees <- as.numeric(selected_parameters$Trees) + 500
selected_eta <- selected_parameters$Eta
selected_gamma <- selected_parameters$Gamma
selected_depth <- selected_parameters$Max_depth
selected_weights <- as.numeric(selected_parameters[c("Weight.loss", "Weight.normal", "Weight.gain")] )

rm(selected_parameters)

# Define the extra parameters
selected_min_child <- 1
selected_seed <- 99

# Get the corresonding rna set
rna_list <- list(
  transcripts_per_million = "tpm",
  scalled_transcripts_per_million = "scld_tpm",
  log_scalled_transcripts_per_million = "log_scld_tpm",
  log_transcripts_per_million = "log_tpm",
  expected_counts = "exp",
  scalled_expected_counts = "scld",
  log_expected_counts = "log",
  log_scalled_expected_counts = "log_scld"
)

rna_selection_name <- rna_list[[selected_rna_set]]

rna_data_path <-  "/hpc/shared/prekovic/dhaynessimmons/data/mRNA/gbm_input/Train/train_"

rna_set <- read.csv(
  paste0(
    rna_data_path,
    rna_selection_name,
    "_soi.csv"
  ),
  row.names = 1
)

# Arm Level Aneuploidies
# Load the data
chr_cnv <- read_tsv(
  "/hpc/shared/prekovic/dhaynessimmons/data/CIN/PANCAN_ArmCallsAndAneuploidyScore_092817.txt"
) %>%
  replace(is.na(.), 0) %>%
  select(-c("Type", "Aneuploidy Score")) %>%
  mutate(Sample = str_replace_all(Sample, "-", "\\.")) %>%
  column_to_rownames("Sample") %>%
  mutate_all(~ replace(., . == 1, 2)) %>%
  mutate_all(~ replace(., . == 0, 1)) %>%
  mutate_all(~ replace(., . == -1, 0))

cat("\n\n CNV size: \n", dim(chr_cnv))

# MODELLING

aneu_cat_metrics_df <- data.frame(
  RNA_Set = character(),
  Trees = numeric(),
  Feature = character(),
  Depth = numeric(),
  Child_weight = numeric(),
  Eta = numeric(),
  Gamma = numeric(),
  Weight_loss = numeric(),
  Weight_norm = numeric(),
  Weight_gain = numeric(),
  Trained_mlogloss = numeric(),
  Test_mlogloss = numeric(),
  Seed = numeric()
)

# Determine the class weights for the target feature

cat("\n selected weights: ")
print(selected_weights)
cat("\n\n")

# Function to map factor levels to weights
feature_digit_function <- function(factors) {
  sapply(factors, function(x) selected_weights[as.numeric(x)])
}

full_df <- merge(rna_set,
  chr_cnv,
  by = "row.names"
)
rm(rna_set)
cat("\n\n full_df: \n")
print(head(full_df[, 1:5]))
cat("\n\n")

y <- as.integer(full_df[[selected_feature]])
X <- full_df %>% select(-c("Row.names", colnames(chr_cnv)))
cat("\n\n Predictors: \n")

# print(head(X[, 1:5]))

train_y_factor <- factor(y, levels = c(0, 1, 2))
weights <- as.numeric(feature_digit_function(train_y_factor))
rm(train_y_factor)

xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y, weight = weights)
rm(weights)
rm(X)
rm(y)

grid <- expand.grid(
  eta = seq(selected_eta - 0.1, selected_eta + 0.1, 0.1),
  depth = seq(selected_depth, selected_depth + 2, 1)
)

set.seed(selected_seed)

for (j in 1:nrow(grid)) {

  for (param in names(grid)) {
    assign(paste0("selected_", param), grid[j, param])
  }

  cat(paste0(
    "\t\t eta: ", selected_eta,
    "\t\t gamma: ", selected_gamma,
    "\t\t depth: ", selected_depth,
    "\t\t trees: ", selected_trees,
    "\t\t child_weight: ", selected_min_child,
    "\n"))

  m_xgb_untuned <- xgb.cv(
    data = xgb_data,
    nrounds = selected_trees,
    objective = "multi:softmax",
    eval_metric = "mlogloss",
    early_stopping_rounds = 100,
    nfold = 5,
    max_depth = selected_depth,
    min_child_weight = selected_min_child,
    eta = selected_eta,
    gamma = selected_gamma,
    num_class = 3,
    stratified = TRUE,
    print_every_n = 25
  )

  best_iteration <- 0

  # First, check if best_iteration is valid
  if (is.null(
    m_xgb_untuned$best_iteration
  ) ||
    m_xgb_untuned$best_iteration < 1) {
    cat(paste0(
      "Warning: No valid best_iteration found.",
      " Using last iteration values instead.\n"
    ))
    # Use the last iteration if best_iteration is not valid
    best_iteration <- nrow(m_xgb_untuned$evaluation_log)
  } else {
    # Ensure that the best_iteration does not exceed the number of rows logged
    if (m_xgb_untuned$best_iteration > nrow(m_xgb_untuned$evaluation_log)) {
      cat(paste0(
        "Warning: best_iteration exceeds the number of rows in evaluation_log.",
        " Adjusting to maximum available.\n"
      ))
      best_iteration <- nrow(m_xgb_untuned$evaluation_log)
    } else {
      best_iteration <- m_xgb_untuned$best_iteration
    }


    best_mlogloss_train <- if (best_iteration > 0) {
      m_xgb_untuned$evaluation_log$train_mlogloss_mean[best_iteration]
    } else {
      NA # Or appropriate default/error value
    }

    best_mlogloss_test <- if (best_iteration > 0) {
      m_xgb_untuned$evaluation_log$test_mlogloss_mean[best_iteration]
    } else {
      NA # Or appropriate default/error value
    }

    cat(paste0(
      "The best iteration occurs with tree #: ",
      best_iteration, "\n\n"
    ))

    aneu_cat_metrics_df <- rbind(aneu_cat_metrics_df, data.frame(
      RNA_Set = selected_rna_set,
      Trees = best_iteration,
      Feature = selected_feature,
      Depth = selected_depth,
      Child_weight = selected_min_child,
      Eta = selected_eta,
      Gamma = selected_gamma,
      Weight_loss = selected_weights[1],
      Weight_norm = selected_weights[2],
      Weight_gain = selected_weights[3],
      Trained_mlogloss = best_mlogloss_train,
      Test_mlogloss = best_mlogloss_test,
      Seed = selected_seed
    ))
  }
}

datetime <- Sys.time() %>%
  str_replace_all(" ", "_") %>%
  str_replace_all(":", "_") %>%
  str_replace_all("\\.", "_")

name <- paste0(
  "/hpc/shared/prekovic/dhaynessimmons/data/model_output/categorical/Cat_xgb_metrics_params_",
  selected_feature, "_",
  index, "_",
  datetime, ".csv"
) %>%
  str_replace_all(" ", "_") %>%
  str_replace_all(":", "_")

write.csv(
  aneu_cat_metrics_df,
  file = name,
  row.names = FALSE
)

cat(paste0("\n Completed processing for index: ", index, "\n"))
