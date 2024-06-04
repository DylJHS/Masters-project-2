library(dplyr)

path <- "/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project/Data/RNA_Data/Model_Input/"
setwd(path)



rna_list <- c(
  "log_scld_tpm_soi.csv",
  "scld_tpm_soi.csv",
  "log_tpm_soi.csv",
  "tpm_soi.csv",
  "exp_soi.csv",
  "log_exp_soi.csv",
  "log_scld_exp_soi.csv",
  "scld_exp_soi.csv"
)

for (i in rna_list) {
  set.seed(99)

  set <- read.csv(
    paste0(path,"Full/", i)
  )

  # Get the total number of rows
  total_rows <- nrow(set)

  # Calculate 70% of the total rows
  num_rows_to_select <- round(0.25 * total_rows)

  # Get a random sample of row indices
  random_indices <- sample(1:total_rows, num_rows_to_select, replace = FALSE)

  train_set <- set %>%
    filter(!row_number() %in% random_indices)

  test_set <- set %>%
    filter(row_number() %in% random_indices)

  test_name <- paste0(path, "test", "_", i)
  train_name <- paste0(path, "train", "_", i)

  # write.csv(
  #   test_set,
  #   test_name,
  #   row.names = FALSE
  # )

  # write.csv(
  #   train_set,
  #   train_name,
  #   row.names = FALSE
  # )
  
}
