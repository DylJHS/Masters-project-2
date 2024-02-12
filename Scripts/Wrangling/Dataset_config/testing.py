import random
import pandas as pd
import time




pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

# Reading the full RNA TPM dataset
# TPM_Samp_df = pd.read_csv("../../../Data/RNA_Data/TCGA_mRNA_TPM.csv")




# GTEx_meta_df = pd.read_csv("../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_Subject_meta.csv")

# # Reading the second metadata dataset that contains the rest of the metadata for the subjects
# GTEx_Annotations_df = pd.read_csv("../../../Data/RNA_Data/GTEx_RNA/GTEx_Annotations.csv")


# Full_tpm_df = pd.read_csv('../../../Data/RNA_Data/TCGA_TPM/TCGA_mRNA_TPM.csv', nrows = 400)
# GTEx_tpm_df = pd.read_csv('../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_TPM/GTEx_RNA_TPM.csv', nrows = 400)

start_time = time.time()

used = set()
t = ['turlte','tle','js','kwj','js',1,2,5,3,5,2,5,2,5,2,5]
used.update(t)

u = []

print(all(elem in t for elem in u))

print(len(used))

end_time = time.time()
duration = end_time - start_time
print("Duration:", duration, "seconds")