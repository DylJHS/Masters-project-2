library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)
library(caret)
library(caTools)

tuning_samples <- 8000
gamma <- 0
depth <- 1

args <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1]) # This is the SLURM_ARRAY_TASK_ID

# RNA SOI SETS
# Expected Counts
exp_set <- read.csv("exp_soi.csv", row.names = 1)[1:tuning_samples,]
cat("\n Exp Count df: \n")
print(head(exp_set[, 1:10]))
scld_exp_set <- read.csv("scld_exp_soi.csv", row.names = 1)[1:tuning_samples,]
log_exp <- read.csv("log_exp_soi.csv", row.names = 1)[1:tuning_samples,]
log_scld_exp <- read.csv("log_scld_exp_soi.csv", row.names = 1)[1:tuning_samples,]

# Transcripts Per Million
tpm_set <- read.csv("tpm_soi.csv", row.names = 1)[1:tuning_samples,]
cat("\n\n TPM df: \n")
print(head(tpm_set[, 1:10]))
scld_tpm_set <- read.csv("scld_tpm_soi.csv", row.names = 1)[1:tuning_samples,]
log_tpm <- read.csv("log_tpm_soi.csv", row.names = 1)[1:tuning_samples,]
log_scld_tpm <- read.csv("log_scld_tpm_soi.csv", row.names = 1)[1:tuning_samples,]

# Pericentromeric CNVs
peri_cnv <- read.csv("../../../../data/CIN/TCGA_pericentro_cnv_hpc.csv")
cat("\n pericentromeric data: \n")
print(head(peri_cnv[, 1:10]))

aneu_feature_list <- colnames(peri_cnv)[-1]

feature <- aneu_feature_list[[index]]
cat(paste0("\n\n feature:", feature, "\n"))

# MODELLING

aneu_reg_metrics_df <- data.frame(
  RNA_Set = character(),
  Feature = character(),
  Depth = numeric(),
  Learning_Rate = numeric(),
  Gamma = numeric(),
  RMSE = numeric()
)

rna_list <- list(
  # transcripts_per_million = tpm_set, # Seems to be the best performing
  scalled_transcripts_per_million = scld_tpm_set, # not too useful (scaled)
  log_scalled_transcripts_per_million = log_scld_tpm,
  # log_transcripts_per_million = log_tpm,
  expected_counts = exp_set,
  scalled_expected_counts = scld_exp_set,
  # log_expected_counts = log_exp,
  log_scalled_expected_counts = log_scld_exp
)


for (i in 1:length(rna_list)) {
  rna <- rna_list[[i]]
  name <- names(rna_list)[i]
  cat(paste0("\t", name, "\n"))
  print(head(rna[, 1:10]))

  full_df <- merge(rna, peri_cnv, by.x = "row.names", by.y = "sampleID")
  y <- as.numeric(full_df[[feature]])
  X <- full_df %>% select(-c("Row.names",all_of(aneu_feature_list))) %>%
     mutate(across(everything(), as.numeric))

  xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y)

  grid <- expand.grid(
    eta = seq(0.01, 0.2, 0.03)
  )

  for (j in 1:nrow(grid)) { # nolint
    cat(paste0(
      "\t\t ETA: ", grid$eta[j], 
      "\n")
      )

    m_xgb_untuned <- xgb.cv(
      data = xgb_data,
      nrounds = 5000,
      objective = "reg:squarederror",
      eval_metric = "rmse",
      early_stopping_rounds = 100,
      nfold = 3,
      max_depth = depth,
      eta = grid$eta[j],
      gamma = gamma,
      verbose = 0
    )

    best_rmse <- m_xgb_untuned$evaluation_log$test_rmse_mean[
      m_xgb_untuned$best_iteration
    ]


    aneu_reg_metrics_df <- rbind(aneu_reg_metrics_df, data.frame(
      Feature = feature,
      RNA_Set = name,
      Depth = depth,
      Learning_Rate = grid$eta[j],
      Gamma = gamma,
      RMSE = best_rmse
    ))
  }
}
write.csv(
  aneu_reg_metrics_df,
  paste0(
    "output/regression/aneu_reg_peri_xgb_metrics_params_", feature, "_", Sys.Date(), ".csv"
  )
)
cat("\n Completed processing for index: ", index, "\n")
