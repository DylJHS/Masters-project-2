library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)
library(caret)
library(caTools)

setwd("/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project")

rna_data_path <- "Data/RNA_Data/Model_Input/Train/train_"
index <- 3
depth <- 5
min_child <- 1
lr <- 0.03

# Default parameters
# trees <- 50
depth <- 5
min_child <- 1
lr <- 0.3

# RNA SOI SETS
# Expected Counts
exp_set <- read.csv(
  paste0(
    rna_data_path,
    "exp_soi.csv"
  ),
  row.names = 1
)

cat("\n Exp Count df: \n")
print(head(exp_set[, 1:10]))

scld_exp_set <- read.csv(
  paste0(
    rna_data_path,
    "scld_exp_soi.csv"
  ),
  row.names = 1
)

log_exp <- read.csv(
  paste0(
    rna_data_path,
    "log_exp_soi.csv"
  ),
  row.names = 1
)

log_scld_exp <- read.csv(
  paste0(
    rna_data_path,
    "log_scld_exp_soi.csv"
  ),
  row.names = 1
)

# Transcripts Per Million
tpm_set <- read.csv(
  paste0(
    rna_data_path,
    "tpm_soi.csv"
  ),
  row.names = 1
)

cat("\n\n TPM df: \n")
print(head(tpm_set[, 1:10]))

scld_tpm_set <- read.csv(
  paste0(
    rna_data_path,
    "scld_tpm_soi.csv"
  ),
  row.names = 1
)

log_tpm <- read.csv(
  paste0(
    rna_data_path,
    "log_tpm_soi.csv"
  ),
  row.names = 1
)

log_scld_tpm <- read.csv(
  paste0(
    rna_data_path,
    "log_scld_tpm_soi.csv"
  ),
  row.names = 1
)

# HRD scores
ori_hrd <- read_tsv("Data/CIN_Features/TCGA.HRD_withSampleID.txt")

# Pericentromeric CNVs
peri_cnv <- read.csv("/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project/Data/CIN_Features/CNV_Data/lim_alpha_incl_TCGA_pericentro.csv") %>%
  mutate_all(~ replace(., is.na(.), 0)) %>%
  mutate(sampleID = gsub("-", ".", sampleID))

cat("\n\n pericentromeric data: \n")
print(head(peri_cnv[, 1:5]))

cat("\n\n All dfs loaded \n")

t_hrd <- as.data.frame(t(ori_hrd))
first_hrd <- t_hrd
colnames(first_hrd) <- t_hrd[1, ]
hrd <- as.data.frame(first_hrd[-1, ]) %>%
  mutate_all(as.numeric) %>%
  rename(loh_hrd = "hrd-loh") %>%
  mutate(new = str_replace_all(rownames(.), "-", "\\.")) %>%
  select(-"HRD")

rownames(hrd) <- hrd$new
hrd <- hrd %>%
  select(-new)

cat("\n\n hrd data: \n")
print(head(hrd))

rm(t_hrd)
rm(first_hrd)

full_cin <- merge(
  hrd,
  peri_cnv,
  by.x = "row.names",
  by.y = "sampleID"
) %>%
  mutate(Row.names = str_replace_all(Row.names, "-", ".")) %>%
  column_to_rownames("Row.names")

cat("\n\n All feature names: ", colnames(full_cin), "\n")

cat("\n\n full cin data: \n")
print(head(full_cin[, 1:5]))

rm(hrd)
rm(peri_cnv)

aneu_reg_feature_list <- colnames(full_cin)

# MODELLING

aneu_reg_metrics_df <- data.frame(
  RNA_Set = character(),
  Trees = numeric(),
  Feature = character(),
  Depth = numeric(),
  Child_weight = numeric(),
  Learning_Rate = numeric(),
  Gamma = numeric(),
  Trained_RMSE = numeric(),
  Test_RMSE = numeric()
)

rna_list <- list(
  transcripts_per_million = tpm_set,
  scalled_transcripts_per_million = scld_tpm_set, # not too useful (scalled)
  log_scalled_transcripts_per_million = log_scld_tpm,
  # log_transcripts_per_million = log_tpm,
  # expected_counts = exp_set,
  # scalled_expected_counts = scld_exp_set,
  # log_expected_counts = log_exp,
  log_scalled_expected_counts = log_scld_exp
)
rna_names <- names(rna_list)

combinations <- expand.grid(
  feature = aneu_reg_feature_list,
  RNA_Set = rna_names,
  stringsAsFactors = FALSE
)

cat("\n\n All combinations: ")
print(combinations)

total_combinations <- nrow(combinations)
cat("\n\n Number fo total combinations: ", total_combinations)

# Select the specific feature and RNA set based on the SLURM task ID
selected_combination <- combinations[index, ]
selected_feature <- selected_combination$feature
selected_rna_set <- selected_combination$RNA_Set

cat(paste0(
  "\n\n Running model for feature: ",
  selected_feature,
  " and RNA set: ",
  selected_rna_set, "\n"
))

# Now select the data based on these choices
rna_data <- rna_list[[selected_rna_set]]
cat("\n\n RNA data: \n")
print(head(rna_data[, 1:5]))
cat("\n\n")

full_df <- merge(rna_data, full_cin, by = "row.names")
cat("\n\n full_df: \n")
print(head(full_df[, 1:5]))

y <- as.numeric(full_df[[selected_feature]])
X <- full_df %>% select(-c("Row.names", colnames(full_cin)))
cat("\n\n Predicotrs: \n")
print(head(X[, 1:5]))

xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y)

grid <- expand.grid(
  gam = seq(0.0, 0.3, 0.5),
  trees = seq(50, 2000, 500)
)

for (j in 1:nrow(grid)) { # nolint
  cat(paste0(
    "\t\t eta: ", lr,
    "\t\t gamma: ", grid$gam[j],
    "\t\t depth: ", depth,
    "\t\t trees: ", grid$trees[j],
    "\t\t child_weight: ", min_child,
    "\n"
  ))

  m_xgb_untuned <- xgb.cv(
    data = xgb_data,
    nrounds = grid$trees[j],
    objective = "reg:squarederror",
    eval_metric = "rmse",
    early_stopping_rounds = 500,
    nfold = 2,
    max_depth = depth,
    min_child_weight = min_child,
    eta = lr,
    gamma = grid$gam[j],
    verbose = 0
  )

  best_iteration <- 0

  # First, check if best_iteration is valid
  if (is.null(
    m_xgb_untuned$best_iteration) ||
    m_xgb_untuned$best_iteration < 1) {
    cat(paste0("Warning: No valid best_iteration found.",
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
  }

  # Accessing the RMSE values safely
  best_rmse_trained <- if (best_iteration > 0) {
    m_xgb_untuned$evaluation_log$train_rmse_mean[best_iteration]
  } else {
    NA # Or appropriate default/error value
  }

  best_rmse_test <- if (best_iteration > 0) {
    m_xgb_untuned$evaluation_log$test_rmse_mean[best_iteration]
  } else {
    NA # Or appropriate default/error value
  }

  cat(paste0(
    "The best iteration occurs with tree #: ",
    best_iteration, "\n\n"
  ))


  aneu_reg_metrics_df <- rbind(aneu_reg_metrics_df, data.frame(
    RNA_Set = selected_rna_set,
    Trees = grid$trees[j],
    Feature = selected_feature,
    Depth = depth,
    Child_weight = min_child,
    Learning_Rate = lr,
    Gamma = grid$gam[j],
    Trained_RMSE = best_rmse_trained,
    Test_RMSE = best_rmse_test
  ))
}


datetime <- Sys.time() %>%
  str_replace_all(" ", "_") %>%
  str_replace_all(":", "_")

name <- paste0(
  "/hpc/shared/prekovic/dhaynessimmons/data/model_output/regression/Reg_xgb_metrics_params_", 
  selected_feature, "_",
  selected_rna_set, "_",
  datetime, ".csv") %>%
  str_replace_all(" ", "_") %>%
  str_replace_all(":", "_")


cat("\n Completed processing for index: ", index, "\n")
