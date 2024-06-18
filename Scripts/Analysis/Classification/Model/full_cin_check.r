
library(dplyr)
library(readr)
library(mltools)
library(knitr)
library(factoextra)
library(psych)
library(tidyverse)

print(Sys.time)
# # HRD scores
# ori_hrd <- read_tsv("Data/CIN_Features/TCGA.HRD_withSampleID.txt")

# # Arm level aneuploidies
# ori_arm_cnv <- read_tsv("Data/CIN_Features/CNV_Data/PANCAN_ArmCallsAndAneuploidyScore_092817.txt") # nolint

# t_hrd <- as.data.frame(t(ori_hrd))
# first_hrd <- t_hrd
# colnames(first_hrd) <- t_hrd[1, ]
# hrd <- as.data.frame(first_hrd[-1, ]) %>%
#   mutate_all(as.numeric) %>%
#   rename(loh_hrd = "hrd-loh")

# # ARM-LEVEL ANEUPLOIDIES
# # Replace the NAs with 0
# cvn_arm <- ori_arm_cnv %>% replace(is.na(.), 0)


# # Re-organise the columns
# cnvs_arm <- cvn_arm %>%
#   column_to_rownames("Sample") %>%
#   dplyr::select(-"Type") %>%
#   mutate_all(as.numeric)


# full_cin <- merge(
#   hrd,
#   cnvs_arm,
#   by = "row.names"
# ) %>%
#   mutate(Row.names = str_replace_all(Row.names, "-", ".")) %>%
#   column_to_rownames("Row.names")

# print(head(full_cin[,1:10]))