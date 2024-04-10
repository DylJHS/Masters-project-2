



# Load the libraries
print("Started")

library(conflicted)
library(dplyr)
library(knitr)
library(tidyverse)
library(edgeR)
library(limma)
library(data.table)

print("Loaded packages \n")

# REUSABLE FUNCTIONS

extract_element <- function(strings, index) {
  # Split each string by "." and extract the third element
  element_list <- sapply(stringr::str_split(strings, "\\."), function(x) x[index])
  return(element_list)
}

untransform <- function(x) { # Function to convert the log transformed counts back into original counts
  return(ceiling((2^x) - 1))
}

# Load the datasets

# RNAseq count data
original <- read.csv("../../../data/tcga_gene_expected_count.csv")

# Metadata
tss_meta <- read.csv("../../../data/meta/tissueSourceSite.tsv", sep = "\t")
abbrv_meta <- read.csv("../../../data/meta/bcrBatchCode.tsv", sep = "\t")
sampletype <- read.csv("../../../data/meta/sampleType.tsv", sep = "\t", colClasses = c("character"))

# Gene Metadata
gene_ids <- read.delim("../../../data/meta/TCGA_PanCan_TPM_Gene_Annotations.txt")

print("Loaded datasets \n")


transformed_data <- original
rownames(transformed_data) <- NULL

real_count_data <- transformed_data %>%
  mutate_at(vars(-1), untransform)

# Format the data

count_data_ordered <- real_count_data[,order(colnames(real_count_data))]

# Combine the metadata
meta <- dplyr::left_join(tss_meta %>% 
                           dplyr::select(c("TSS.Code", "Study.Name")) %>% 
                           distinct() %>% 
                           sapply(trimws) %>% 
                           as.data.frame(),
                         abbrv_meta %>%
                           dplyr::select(c("Study.Abbreviation", "Study.Name")) %>% 
                           distinct()%>% 
                           sapply(trimws) %>% 
                           as.data.frame(), 
                         by = "Stud                 y.Name")

# CONVERT THE GENE IDS INTO GENE NAMES

counts_data <- dplyr::right_join(gene_ids %>% 
                    dplyr::select(c("id", "gene")) %>% 
                    sapply(trimws) %>% 
                    as.data.frame(),
                  count_data_ordered,
                  by = c("id" = "sample")) %>% 
  dplyr::select(-"id")

count_data <- counts_data %>%
  mutate(gene = trimws(gene))

print("\n count_data")
print(head(count_data[, 1:10]))


# 2.  Edit the samples (column names) into ids so as to be able to identify the participant and the sample and cancer type

# 2.1. extract the different elements from the sample names

ids <- colnames(count_data)[-1]
participants <- extract_element(ids, 3)
condition <- extract_element(ids, 4)
tissue_type <- extract_element(ids, 2)

# 2.2. Reform the column names using the extracted elements and reset the index as the Genes

column_names <- c("Genes" ,  paste(tissue_type, participants, condition, sep = "."))
colnames(count_data) <- column_names

# Convert data to a data.table for faster processing during the grouping of the duplicate genes
count_data <- data.table::setDT(count_data)

# Combine duplicate genes together using the median of the expression
grouped <- count_data[, lapply(.SD, function(x) if (length(x) > 1) ceiling(median(x, na.rm = TRUE)) else x), by = Genes, .SDcols = -"Genes"]

data <- grouped %>%
  column_to_rownames(var = "Genes")

print("\n grouped")
print(head(grouped[, 1:10]))


# Remove the unwanted sample types from the samples based on the code
codes_to_use <- c("01", "02", "03", "04", "05", "08", "09")
samples_to_use <- data %>%  dplyr::select(ends_with(codes_to_use))

# 2.3. Assign the TSS abbreviation to the column name in place of the TSS code

column_names <- colnames(samples_to_use)
prefixes <- substr(column_names,1,2)
abbrvs <- meta$Study.Abbreviation[match(prefixes, meta$TSS.Code)]
column_names_prefixes <- paste0(abbrvs, substr(column_names, 3, nchar(column_names)))
colnames(samples_to_use) <- column_names_prefixes

# 2.4. Assign the sample type abbreviation to the column name in place of the sample type code

new_column_names <- colnames(samples_to_use)
suffixes <- substr(new_column_names, nchar(new_column_names) - 1, nchar(new_column_names))
s_abbrvs <- sampletype$Short.Letter.Code[match(suffixes, sampletype$Code)]
column_names_suffxes <- paste0(substr(new_column_names, 1, nchar(new_column_names) - 2), s_abbrvs)
colnames(samples_to_use) <- column_names_suffxes

print("\n samples_to_use")
print(head(samples_to_use[, 1:10]))


# 3.1. Create the tissue & sample interactions

# sampled_data <- as.data.frame(colnames(samples_to_use)) %>% 
#   dplyr::group_by(TSS := extract_element(colnames(samples_to_use), 1),
#                    ST := extract_element(colnames(samples_to_use), 3)) %>% 
#   dplyr::select(TSS, ST ) %>% 
#   dplyr::count()

# # 3.2. Find the groups that have less than 10 samples

# less_than10 <- sampled_data %>% 
#   dplyr::filter(freq < 10) 

# 4.1. DGE MODELLING: DEFINE THE MATRIX, GET THE NORMALISATION FACTORS AND SCALE THE DATA

# convert df to matrix 
mapped_matrix <- as.matrix(samples_to_use)

# Calculate the normalisation factor 
d0 <- edgeR::calcNormFactors(mapped_matrix, method = "TMM")

# Scale Gene counts

# Divide the samples by their corresponding normalisation factors
scaled_data <- base::sweep(samples_to_use, 2, d0, "/")
scaled_data <- ceiling(scaled_data) # Round up

# Remove 0 genes and samples
unzero_d <- scaled_data[rowSums(scaled_data[]) > 0, ]
cols_zeros <- which(apply(unzero_d, 2, function(x) all(x == 0)))
if (length(cols_zeros) > 0) {
  reduced_scaled <- unzero_d[, -cols_zeros]
} else {
  reduced_scaled <- unzero_d
}

# 5.  MDS

# Transpose the data 
t_scaled_data <- t(reduced_scaled)
write.csv(t_scaled_data, "Transposed_scaled_&_reduced_df.csv", row.names = TRUE)

print("\n t_scaled_data")
print(head(t_scaled_data[, 1:3]))

sample_dist_matrix <- dist(t_scaled_data)
sample_dist_df <- as.data.frame(as.matrix(sample_dist_matrix))

print("\n sample_dist_matrix")
print(head(sample_dist_df[, 1:10]))

write.csv(sample_dist_df, "sample_dist_matrix.csv", row.names = TRUE)

