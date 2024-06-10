library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)
library(caret)
library(caTools)


args <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1]) # This is the SLURM_ARRAY_TASK_ID

# Default parameters
trees <- 50
depth <- 5

cat("Index: ", index, "\n\n")

rna_data_path <- "/hpc/shared/prekovic/dhaynessimmons/data/mRNA/gbm_input/Train/train_"

# RNA SOI SETS
# Expected Counts
exp_set <- read.csv(
  paste0(
    rna_data_path,
    "exp_soi.csv"
  ),
  row.names = 1
)

cat("\n Exp Count df: \n\n")
cat(sprintf("   %s", head(colnames(exp_set), 10)), "\n", sep = "")
cat(sprintf("Dimensions: %d rows, %d columns\n\n", nrow(exp_set), ncol(exp_set)), "\n")

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

cat("\n\n TPM df: \n\n")
cat(sprintf("   %s", head(colnames(tpm_set), 10)), "\n", sep = "")
cat(sprintf("Dimensions: %d rows, %d columns\n\n", nrow(tpm_set), ncol(tpm_set)), "\n\n")

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
ori_hrd <- read_tsv("/hpc/shared/prekovic/dhaynessimmons/data/CIN/TCGA.HRD_withSampleID.txt")

# Pericentromeric CNVs
peri_cnv <- read.csv(
  "/hpc/shared/prekovic/dhaynessimmons/data/CIN/lim_alpha_incl_TCGA_pericentro_cnv.csv"
) %>%
  mutate_all(~ replace(., is.na(.), 0)) %>%
  mutate(sampleID = gsub("-", ".", sampleID))

cat("\n pericentromeric data: \n")
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
  RMSE = numeric()
)

rna_list <- list(
  transcripts_per_million = tpm_set,
  scalled_transcripts_per_million = scld_tpm_set, # not too useful (scalled)
  log_scalled_transcripts_per_million = log_scld_tpm,
  # log_transcripts_per_million = log_tpm,
  # expected_counts = exp_set
  # scalled_expected_counts = scld_exp_set,
  # log_expected_counts = log_exp,
  log_scalled_expected_counts = log_scld_exp
)
rna_names <- names(rna_list)

combinations <- expand.grid(feature = aneu_reg_feature_list, RNA_Set = rna_names, stringsAsFactors = FALSE)
cat("\n\n All combinations: ")
print(combinations)

total_combinations <- nrow(combinations)
cat("\n\n Number of total combinations: ", total_combinations)

# Select the specific feature and RNA set based on the SLURM task ID
selected_combination <- combinations[index, ]
selected_feature <- selected_combination$feature
selected_rna_set <- selected_combination$RNA_Set

cat(paste0("\n\n Running model for feature: ", selected_feature, " and RNA set: ", selected_rna_set, "\n"))

# Now select the data based on these choices
rna_data <- rna_list[[selected_rna_set]]
cat("\n\n RNA data: \n")
print(head(rna_data[, 1:5]))


full_df <- merge(rna_data, full_cin, by = "row.names")
y <- as.numeric(full_df[[selected_feature]])
X <- full_df %>% select(-c("Row.names", colnames(full_cin)))

xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y)

grid <- expand.grid(
  lr = seq(0.25, 0.29, 0.2),
  gam = seq(0, 0.15, 0.2),
  child_weight = seq(0.3, 2.1, 0.3)
)

for (j in 1:nrow(grid)) { # nolint
  cat(paste0(
    "\t\t eta: ", grid$lr[j],
    "\t\t gamma: ", grid$gam[j],
    "\t\t depth: ", depth,
    "\t\t trees: ", trees,
    "\t\t child_weight: ", grid$child_weight[j],
    "\n"
  ))

  m_xgb_untuned <- xgb.cv(
    data = xgb_data,
    nrounds = 10000,
    objective = "reg:squarederror",
    eval_metric = "rmse",
    early_stopping_rounds = 50,
    nfold = 10,
    max_depth = depth,
    eta = grid$lr[j],
    gamma = grid$gam[j],
    min_child_weight = grid$child_weight[j],
    n_estimators = trees,
    verbose = 0
  )

  best_rmse <- m_xgb_untuned$evaluation_log$test_rmse_mean[
    m_xgb_untuned$best_iteration
  ]

  aneu_reg_metrics_df <- rbind(aneu_reg_metrics_df, data.frame(
    Feature = selected_feature,
    RNA_Set = selected_rna_set,
    Trees = trees,
    Depth = depth,
    Child_weight = grid$child_weight[j],
    Learning_Rate = grid$lr[j],
    Gamma = grid$gam[j],
    RMSE = best_rmse
  ))
}

datetime <- Sys.time() %>%
  str_replace_all("[ :.]", "_")

name <- paste0(
  "/hpc/shared/prekovic/dhaynessimmons/data/model_output/regression/Reg_xgb_metrics_params_", selected_feature, "_", selected_rna_set, "_", datetime, ".csv"
) %>%
  str_replace_all("[ :]", "_")

write.csv(
  aneu_reg_metrics_df,
  file = name,
  row.names = FALSE
)

cat("\n Completed processing for index: ", index, "\n")
