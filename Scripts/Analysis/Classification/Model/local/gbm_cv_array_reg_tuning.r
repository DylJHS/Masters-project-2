library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)
library(caret)
library(caTools)

setwd("/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project")

rna_data_path <- "Data/RNA_Data/Model_Input/Train/train_"


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
peri_cnv <- read.csv("Data/CIN_Features/CNV_Data/TCGA_pericentro_cnv_hpc.csv")
cat("\n\n pericentromeric data: \n")
print(head(peri_cnv[, 1:10]))

cat("\n\n All dfs loaded \n")

t_hrd <- as.data.frame(t(ori_hrd))
first_hrd <- t_hrd
colnames(first_hrd) <- t_hrd[1, ]
hrd <- as.data.frame(first_hrd[-1, ]) %>%
  mutate_all(as.numeric) %>%
  rename(loh_hrd = "hrd-loh") %>%
  mutate(new = str_replace_all(rownames(.), "-", "\\."))

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

cat("\n\n full cin data: \n")
print(head(full_cin))

rm(hrd)
rm(peri_cnv)

aneu_reg_feature_list <- colnames(full_cin)

feature <- aneu_reg_feature_list[[index]]
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
  transcripts_per_million = tpm_set,
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
  y <- as.numeric(full_df[[feature]])
  X <- full_df %>% select(-c("Row.names", colnames(full_cin)))

  xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y)

  grid <- expand.grid(
    lr = seq(0.05, 0.2, 0.05),
    gam = seq(0, 0.3, 0.2),
    depth = seq(1, 3, 1)
  )


  for (j in 1:nrow(grid)) { # nolint
    cat(paste0(
      "\t\t eta: ", grid$lr[j],
      "\t\t gamma: ", grid$gam[j],
      "\t\t depth: ", grid$depth[j],
      "\n"
    ))

    m_xgb_untuned <- xgb.cv(
      data = xgb_data,
      nrounds = 5000,
      objective = "reg:squarederror",
      eval_metric = "rmse",
      early_stopping_rounds = 50,
      nfold = 3,
      max_depth = grid$depth[j],
      eta = grid$lr[j],
      gamma = grid$gam[j],
      verbose = 0
    )

    best_rmse <- m_xgb_untuned$evaluation_log$test_rmse_mean[
      m_xgb_untuned$best_iteration
    ]


    aneu_reg_metrics_df <- rbind(aneu_reg_metrics_df, data.frame(
      Feature = feature,
      RNA_Set = name,
      Depth = grid$depth[j],
      Learning_Rate = grid$lr[j],
      Gamma = grid$gam[j],
      RMSE = best_rmse
    ))
  }
}
write.csv(
  aneu_reg_metrics_df,
  paste0(
    "output/regression/aneu_reg_xgb_metrics_params_", feature, "_", Sys.time(), ".csv"
  ),
  row.names = FALSE
)
cat("\n Completed processing for index: ", index, "\n")
