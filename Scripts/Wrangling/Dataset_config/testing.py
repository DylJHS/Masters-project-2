import random
import pandas as pd



pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

# Reading the full RNA TPM dataset
# TPM_Samp_df = pd.read_csv("../../../Data/RNA_Data/TCGA_mRNA_TPM.csv")




GTEx_meta_df = pd.read_csv("../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_Subject_meta.csv")

# Reading the second metadata dataset that contains the rest of the metadata for the subjects
GTEx_Annotations_df = pd.read_csv("../../../Data/RNA_Data/GTEx_RNA/GTEx_Annotations.csv")


raw_GTEx_tpm_df = pd.read_csv('../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_TPM/Raw_GTEx_RNA_TPM.txt', sep='\t', skiprows=2)


print(raw_GTEx_tpm_df.iloc[:5,:5])





