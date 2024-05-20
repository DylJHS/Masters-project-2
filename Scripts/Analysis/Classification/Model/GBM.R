library(dplyr)
library(readr)
library(mltools)
library(knitr)
library(factoextra)
library(data.table)
library(ggpubr)
library(psych)
library(limma)
library(edgeR)
library(tidyverse)
library(xgboost)
library(caTools)
library(dplyr)
library(caret)

extract_element <- function(strings, index) {
  # Split each string by "." and extract the third element
  element_list <- sapply(strsplit(strings, "\\."), function(x) x[index])
  return(element_list)
}

untransform_exp <- function(x) {#Function to convert the log transformed counts back into original counts
  return(ceiling((2^x)-1))
}

transform_exp <- function(x) {#Function to convert the log transformed counts back into original counts
  return(log2(x + 1))
}

untransform_tpm <- function(x) {#Function to convert the log transformed counts back into original counts
  return(ceiling((2^x)-0.001))
}

transform_tpm <- function(x) {#Function to convert the log transformed counts back into original counts
  return(log2(x + 0.001))
}

# SOI genes
soi <- read.csv("../../../data/mRNA/TCGA_mRNA_TPM_SOI.csv")
soi_genes <-soi[,2]

# Predictor variables

# Expected Counts
ori_exp <- read.csv("../../../data/mRNA/tcga_gene_expected_count.csv")
order_exp <- ori_exp[,order(colnames(ori_exp))]


# TPM counts
ori_tpm <- read.csv("../../../data/mRNA/TCGA_mRNA_TPM_Full.csv")
order_tpm <- ori_tpm[,order(colnames(ori_tpm))] %>% 
  dplyr::select(-"id")



# Number of random rows and columns to select
n_rows <- 30000
row_samples <- sample(nrow(order_exp), n_rows)


order_tpm <- order_tpm[row_samples,0:5000]
order_exp <- order_exp[row_samples,0:5000]

# Response data

# HRD scores
ori_hrd <- read_tsv("../../../data/CIN/TCGA.HRD_withSampleID.txt")

# Arm level aneuploidies
ori_arm_cnv <- read_tsv("../../../data/CIN/PANCAN_ArmCallsAndAneuploidyScore_092817.txt")


# Metadata

# Sample meta
tss_meta <- read.csv("../../../data/meta/tissueSourceSite.tsv", sep = "\t")
abbrv_meta <- read.csv("../../../data/meta/bcrBatchCode.tsv", sep = "\t")

# Gene Metadata
gene_ids <- read.delim("../../../data/meta/TCGA_PanCan_TPM_Gene_Annotations.txt")

meta <- left_join(tss_meta %>% 
                    dplyr::select(c("TSS.Code", "Study.Name")) %>% 
                    distinct() %>% 
                    sapply(trimws) %>% 
                    as.data.frame(),
                  abbrv_meta %>%
                    dplyr::select(c("Study.Abbreviation", "Study.Name")) %>% 
                    distinct()%>% 
                    sapply(trimws) %>% 
                    as.data.frame(), 
                  by = "Study.Name")

trans_exp <- order_exp # Log transformed Expected Counts
rownames(trans_exp) <- NULL
trans_tpm <- order_tpm # Log transformed TPM Counts
rownames(trans_tpm) <- NULL

count_exp <- trans_exp %>%
  mutate_at(vars(-1), untransform_exp)
count_tpm <- trans_tpm %>%
  mutate_at(vars(-1), untransform_tpm)

# Remove the non-cancerous sample types from the set
codes_to_use <- c("01","02","03","04","05","08","09")

exp_samples_to_use <- count_exp %>%  dplyr::select(c("sample", ends_with(codes_to_use)))
tpm_samples_to_use <- count_tpm %>%  dplyr::select(c("Gene", ends_with(codes_to_use)))

# Convert the Gene Ids into names
counts_exp <- right_join(gene_ids %>% 
                           dplyr::select(c("id", "gene")) %>% 
                           sapply(trimws) %>% 
                           as.data.frame(),
                         exp_samples_to_use,
                         by = c("id" = "sample")) %>% 
  dplyr::select(-"id")

counts_exp <- counts_exp %>%
  mutate(gene = trimws(gene))

# Convert data to a data.table for faster processing during the grouping of the duplicate genes
setDT(counts_exp)  
setDT(tpm_samples_to_use)  

# Combine duplicate genes together using the median of the expression
grouped_exp <- counts_exp[, lapply(.SD, function(x) if (length(x) > 1) median(x, na.rm = TRUE) else x), by = gene, .SDcols = -"gene"]
grouped_tpm <- tpm_samples_to_use[, lapply(.SD, function(x) if (length(x) > 1) median(x, na.rm = TRUE) else x), by = Gene, .SDcols = -"Gene"]

# Reset the row names to the Gene names
groups_exp <- as.data.frame(grouped_exp)
groups_tpm <- as.data.frame(grouped_tpm)

exp_data <- distinct(groups_exp) %>% 
  column_to_rownames(var = "gene")
tpm_data <- distinct(groups_tpm) %>% 
  column_to_rownames(var = "Gene")

# Remove samples (cols) with 0 values throughout
cols_zeros_tpm <- which(apply(tpm_data, 2, function(x) all(x == 0)))
cols_zeros_exp <- which(apply(exp_data, 2, function(x) all(x == 0)))

if (length(cols_zeros_tpm) > 0) {
  data_complete_tpm <- tpm_data[, -cols_zeros_tpm]
} else {
  data_complete_tpm <- tpm_data
}

if (length(cols_zeros_exp) > 0) {
  data_complete_exp <- exp_data[, -cols_zeros_exp]
} else {
  data_complete_exp <- exp_data
}

# Identify genes with 0 expression across all samples
zero_genes_exp <- rowSums(data_complete_exp == 0) == ncol(data_complete_exp)

# Remove genes with 0 expression across all samples
filt_dt_exp <- data_complete_exp[!zero_genes_exp, ]

# Turn into a DGE object
matrix_exp <- filt_dt_exp %>% as.matrix()
d_exp <- DGEList(matrix_exp)

# Calculate the normalisation factor 
Normfact_exp <- calcNormFactors(d_exp, method = "TMM")
Normfactors <- as.data.frame(Normfact_exp$samples) %>% select("norm.factors")

# match the column names from the normalisation factor df with that of the other dfs in order to divide the count values by the library size factors for the correct samples 
matching_cols_tpm <- colnames(tpm_data)[match(rownames(Normfactors), colnames(tpm_data))] 

matching_cols_tpm <- matching_cols_tpm[!is.na(matching_cols_tpm)]

scld_cnts_tpm <- round(sweep(tpm_data[, matching_cols_tpm], 2, Normfactors$norm.factors, "/"),2)

matching_cols_exp <- colnames(exp_data)[match(rownames(Normfactors), colnames(exp_data))]
matching_cols_exp <- matching_cols_exp[!is.na(matching_cols_exp)]

scld_cnts_exp <- round(sweep(exp_data[, matching_cols_exp], 2, Normfactors$norm.factors, "/"),2)

# Expected Counts
exp_set <- data_complete_exp %>% # Raw expected counts
  filter(rownames(.) %in% soi_genes) %>%
  t() %>% 
  as.data.frame()

scld_exp_set <- scld_cnts_exp %>% # Library size normalised Expected counts
  filter(rownames(.) %in% soi_genes) %>% 
  t() %>% 
  as.data.frame()

log_exp <- exp_set %>% # Log transformed expected counts
  mutate_at(vars(everything()), transform_exp) 

log_scld_exp <- scld_exp_set %>% # Log transformed Library size normalised expected counts
  mutate_at(vars(everything()), transform_exp) 

# Transcript per million (TPM)
tpm_set <- data_complete_tpm %>% # Raw TPM counts 
  filter(rownames(.) %in% soi_genes) %>% 
  t() %>% 
  as.data.frame()

scld_tpm_set <- scld_cnts_tpm %>% # Library size normalised TPM counts
  filter(rownames(.) %in% soi_genes) %>% 
  t() %>% 
  as.data.frame()

log_tpm <- tpm_set %>% # Log transformed TPM counts
  mutate_at(vars(everything()), transform_tpm) 

log_scld_tpm <- scld_tpm_set %>% # Log transformed Library size normalised TPM counts
  mutate_at(vars(everything()), transform_tpm) 

t_hrd <- as.data.frame(t(ori_hrd)) 
first_hrd <- t_hrd
colnames(first_hrd) <- t_hrd[1, ]
hrd <- as.data.frame(first_hrd[-1,]) %>% 
  mutate_all(as.numeric) %>% 
  rename(loh_hrd = "hrd-loh")


# ARM-LEVEL ANEUPLOIDIES
# Replace the NAs with 0
cvn_arm <- ori_arm_cnv %>% replace(is.na(.), 0)

# Re-organise the columns 
cnvs_arm <- cvn_arm %>%
  column_to_rownames("Sample") %>% 
  dplyr::select(-"Type") %>% 
  mutate_all(as.numeric)

full_cin <- merge(
  hrd,
  cnvs_arm,
  by = "row.names"
) %>% 
  mutate(Row.names = str_replace_all(Row.names, "-", ".")) %>%
  column_to_rownames("Row.names")

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

# aneu_feature_list <- colnames(full_cin[1,6:length(full_cin)])
aneu_feature_list <- colnames(full_cin[1,6:12])

aneu_cat_metrics_df <- data.frame(
  RNA_Set = character(),
  Feature = character(),
  Depth = numeric(),
  Learning_Rate = numeric(),
  Gamma = numeric(),
  Accuracy = numeric()
)

lr= 0.04
set.seed(33)

for (feature in aneu_feature_list){
  cat(paste0("\n", feature, ":"))
  
  for (i in 1:length(rna_list)){
    
    rna <- rna_list[[i]]
    name <- names(rna_list)[i]
    
    cat(paste0("\n\t", name, ":"))
    
    full_df <- merge(rna, full_cin, by = "row.names")
    
    sample_split <- sample.split(Y = full_df[colnames(full_cin)], SplitRatio = 0.7)
    train_set <- subset(x = full_df, sample_split == TRUE)
    test_set <- subset(x = full_df, sample_split == FALSE)
    
    
    y_train <- as.integer(train_set[[feature]])
    X_train <- train_set %>% select(-append("Row.names",colnames(full_cin)))
    
    y_test <- as.integer(test_set[[feature]])
    X_test <- test_set %>% select(-append("Row.names",colnames(full_cin)))
    
    y_train[y_train == -1] <- 0
    y_train[y_train == 1] <- 2
    y_train[y_train == 0] <- 1
    
    y_test[y_test == -1] <- 0
    y_test[y_test == 1] <- 2
    y_test[y_test == 0] <- 1
    
    
    xgb_train <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
    xgb_test <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)
    
    watchlist <- list(train = xgb_train, test = xgb_test)
    
    for (depth in seq(1,14,2)){
      cat(paste0("\n\t\t Depth: ", depth, "\n"))
      for (gam in seq(0.5,3.5,1.5)) {
        cat(paste0("\t\t\t Gamma: ", gam))
        
        xgb_params <- list(
          booster = "gbtree",
          num_class=3,
          max_depth = depth,
          gamma = gam,
          subsample = 0.75,
          colsample_bytree = 1,
          objective = "multi:softmax",
          eval_metric = "mlogloss", 
          eta = 0.04 # Seems like the best eta is the closest one to 0
        )
        
        xgb_model <- xgb.train(
          params = xgb_params,
          data = xgb_train,
          watchlist = watchlist, 
          early_stopping_rounds = 50,
          nrounds = 6000,
          verbose = 0
        )
        
        # Best loss and iteration number
        best_loss <- min(xgb_model[["evaluation_log"]][["test_mlogloss"]])
        
        xgb_cat_preds <- predict(xgb_model, as.matrix(X_test), reshape = TRUE)
        
        xgb_cat_preds <- as.data.frame(xgb_cat_preds)
        xgb_cat_preds <- xgb_cat_preds %>% 
          rename(pred = "xgb_cat_preds") %>% 
          mutate(actual = y_test,
                 correct = ifelse(pred == actual, 1, 0)
          )
        
        aneu_cat_metrics_df <- rbind(aneu_cat_metrics_df,
                                     data.frame(
                                       RNA_Set = name,
                                       Feature = feature,
                                       Depth = depth,
                                       Learning_Rate = lr,
                                       Gamma = gam,
                                       Accuracy = 
                                         sum(xgb_cat_preds$correct)/
                                         length(xgb_cat_preds$correct)                                  )
        )
      }
    }
  }
  write.csv(aneu_cat_metrics_df, paste0("../../../data/model_output/aneu_cat_xgb_metrics_params_",feature,"_",Sys.Date(), ".csv"))
  
}

reg_feature_list <- colnames(full_cin[1,6:12])

set.seed(35)

lr = 0.04


for (feature in reg_feature_list){
  cat(paste0("\n", feature, ":"))
  
  aneu_reg_metrics_df <- data.frame(
    RNA_Set = character(),
    Feature = character(),
    Depth = numeric(),
    Learning_Rate = numeric(),
    Gamma = numeric(),
    Accuracy = numeric()
  )
  
  for (i in 1:length(rna_list)){
    
    rna <- rna_list[[i]]
    name <- names(rna_list)[i]
    cat(paste0("\n\t", name, ": \n"))
    
    full_df <- merge(rna, full_cin, by = "row.names")
    
    sample_split <- sample.split(Y = full_df[colnames(full_cin)], SplitRatio = 0.7)
    train_set <- subset(x = full_df, sample_split == TRUE)
    test_set <- subset(x = full_df, sample_split == FALSE)
    
    y_train <- as.integer(train_set[[feature]])
    X_train <- train_set %>% select(-append("Row.names",colnames(full_cin)))
    
    y_test <- as.integer(test_set[[feature]])
    X_test <- test_set %>% select(-append("Row.names",colnames(full_cin)))
    
    xgb_train <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
    xgb_test <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)
    
    watchlist <- list(train = xgb_train, test = xgb_test)
    
    for (depth in seq(1,14,2)){
      cat(paste0("\n\t\t Depth: ", depth, "\n"))
      # cat(paste0("\t\t Learning Rate: ", 0.04, "\n"))
      for (gam in seq(0.5,3.5,1.5)){
        cat(paste0("\t\t\t Gamma: ", gam))
        xgb_params <- list(
          booster = "gbtree",
          max_depth = depth,
          gamma = gam,
          subsample = 0.75,
          colsample_bytree = 1,
          objective = "reg:squarederror",
          eval_metric = "rmse", 
          eta = lr # Seems like the best eta is the closest one to 0
        )
        
        xgb_model <- xgb.train(
          params = xgb_params,
          data = xgb_train,
          watchlist = watchlist, 
          early_stopping_rounds = 50,
          nrounds = 6000,
          verbose = 0
        )
        
        # xgb_model
        
        # Best RMSE and iteration number
        best_rmse <- min(xgb_model[["evaluation_log"]][["test_rmse"]])
        
        xgb_preds <- predict(xgb_model, as.matrix(X_test), reshape = TRUE)
        xgb_preds <- as.data.frame(xgb_preds)
        xgb_preds <- xgb_preds %>% 
          rename(base_pred = "xgb_preds") %>% 
          mutate(mod_pred = round(base_pred),
                 actual = y_test,
                 correct = ifelse(mod_pred == actual, 1, 0)
          )
        
        aneu_reg_metrics_df <- rbind(aneu_reg_metrics_df,
                                     data.frame(
                                       RNA_Set = name,
                                       Feature = feature,
                                       Depth = depth,
                                       Learning_Rate = lr,
                                       Gamma = gam,
                                       Accuracy = 
                                         sum(xgb_preds$correct)/
                                         length(xgb_preds$correct)
                                       
                                     )
        )
      }
      
    }
    
  }
  write.csv(aneu_reg_metrics_df, paste0("../../../data/model_output/aneu_reg_xgb_metrics_params_",feature,"_",Sys.Date(), ".csv"))
}


