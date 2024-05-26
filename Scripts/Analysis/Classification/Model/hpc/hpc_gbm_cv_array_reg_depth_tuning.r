library(dplyr)
library(readr)
library(mltools)
library(knitr)
library(factoextra)
library(psych)
library(tidyverse)
library(xgboost)
library(caret)
library(caTools)

tuning_samples <- 8000
gamma <- 0
eta <- 0.3

args <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1]) # This is the SLURM_ARRAY_TASK_ID

# RNA SOI SETS
# Expected Counts
exp_set <- read.csv("exp_soi.csv", row.names = 1)[1:8000,]
cat("\n Exp Count df: \n")
print(head(exp_set[, 1:10]))
scld_exp_set <- read.csv("scld_exp_soi.csv", row.names = 1)[1:8000,]
log_exp <- read.csv("log_exp_soi.csv", row.names = 1)[1:8000,]
log_scld_exp <- read.csv("log_scld_exp_soi.csv", row.names = 1)[1:8000,]

# Transcripts Per Million
tpm_set <- read.csv("tpm_soi.csv", row.names = 1)[1:8000,]
cat("\n\n TPM df: \n")
print(head(tpm_set[, 1:10]))
scld_tpm_set <- read.csv("scld_tpm_soi.csv", row.names = 1)[1:8000,]
log_tpm <- read.csv("log_tpm_soi.csv", row.names = 1)[1:8000,]
log_scld_tpm <- read.csv("log_scld_tpm_soi.csv", row.names = 1)[1:8000,]


# HRD scores
ori_hrd <- read_tsv("../../../../data/CIN/TCGA.HRD_withSampleID.txt")

# Arm level aneuploidies
ori_arm_cnv <- read_tsv("../../../../data/CIN/PANCAN_ArmCallsAndAneuploidyScore_092817.txt")
cat("\n All dfs loaded \n")

t_hrd <- as.data.frame(t(ori_hrd))
first_hrd <- t_hrd
colnames(first_hrd) <- t_hrd[1, ]
hrd <- as.data.frame(first_hrd[-1, ]) %>%
  mutate_all(as.numeric) %>%
  rename(loh_hrd = "hrd-loh")
rm(t_hrd)
rm(first_hrd)

# ARM-LEVEL ANEUPLOIDIES
# Replace the NAs with 0
cvn_arm <- ori_arm_cnv %>% replace(is.na(.), 0)
rm(ori_arm_cnv)

# Re-organise the columns
cnvs_arm <- cvn_arm %>%
  column_to_rownames("Sample") %>%
  dplyr::select(-"Type") %>%
  mutate_all(as.numeric)
rm(cvn_arm)

full_cin <- merge(
  hrd,
  cnvs_arm,
  by = "row.names"
) %>%
  mutate(Row.names = str_replace_all(Row.names, "-", ".")) %>%
  column_to_rownames("Row.names")
rm(hrd)
rm(cnvs_arm)

aneu_feature_list <- colnames(full_cin[1, 1:5])

feature <- aneu_feature_list[[index]]
cat(paste0("\n feature:", feature, "\n"))

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
  transcripts_per_million = tpm_set, # Seems to be the best performing
  scalled_transcripts_per_million = scld_tpm_set, # not too useful (scalled)
  log_scalled_transcripts_per_million = log_scld_tpm,
  log_transcripts_per_million = log_tpm,
  expected_counts = exp_set,
  scalled_expected_counts = scld_exp_set,
  log_expected_counts = log_exp,
  log_scalled_expected_counts = log_scld_exp
)


for (i in 1:length(rna_list)) {
  rna <- rna_list[[i]]
  name <- names(rna_list)[i]
  cat(paste0("\t", name, "\n"))

  full_df <- merge(rna, full_cin, by = "row.names")
  y <- as.integer(full_df[[feature]])
  X <- full_df %>% select(-c("Row.names",colnames(full_cin)))

  xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y)

  grid <- expand.grid(
    max_depth = seq(1, 10, 2)
  )

  for (j in 1:nrow(grid)) { # nolint
    cat(paste0(
      "\t\t Depth: ", grid$max_depth[j], 
      "\n")
      )

    m_xgb_untuned <- xgb.cv(
      data = xgb_data,
      nrounds = 5000,
      objective = "reg:squarederror",
      eval_metric = "rmse",
      early_stopping_rounds = 50,
      nfold = 3,
      max_depth = grid$max_depth[j],
      eta = eta,
      gamma = gamma,
      num_class = 3,
      verbose = 0
    )

    best_rmse <- m_xgb_untuned$evaluation_log$test_rmse_mean[
      m_xgb_untuned$best_iteration
    ]


    aneu_reg_metrics_df <- rbind(aneu_reg_metrics_df, data.frame(
      Feature = feature,
      RNA_Set = name,
      Depth = grid$max_depth[j],
      Learning_Rate = eta,
      Gamma = gamma,
      RMSE = best_rmse
    ))
  }
}
write.csv(
  aneu_reg_metrics_df,
  paste0(
    "output/regression/aneu_reg_xgb_metrics_params_", feature, "_", Sys.Date(), ".csv"
  )
)
cat("\n Completed processing for index: ", index, "\n")
