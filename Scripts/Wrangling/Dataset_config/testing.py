import random
import pandas as pd
import time

pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

# GTEx_meta_df = pd.read_csv("../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_Subject_meta.csv")
GTEx_tpm_df = pd.read_csv('../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_TPM/GTEx_RNA_TPM.csv', nrows = 500 )
# GTEx_Annotations_df = pd.read_csv("../../../Data/RNA_Data/GTEx_RNA/GTEx_Annotations.csv")

TPM_Samp_df = pd.read_csv("../../../Data/RNA_Data/TCGA_TPM/TCGA_mRNA_TPM.csv", nrows = 400)
# TCGA_Sample_type = pd.read_csv("../../../Data/Other/TCGA_meta/sampleType.tsv", sep = '\t')


n=1
while True:
	Normal_set = pd.read_csv(f'../../../Data/RNA_Data/TCGA_TPM/TCGA_Normal_CTRL_Sets/TCGA_Normal_mRNA_TPM_CTRL_Set{n}.csv')
	n +=1
	print(Normal_set.shape)
Control_set = pd.read_csv('../../../Data/RNA_Data/TCGA_TPM/TCGA_Cancer_CTRL_Sets/TCGA_TPM_RNA_Control_df1.csv')




