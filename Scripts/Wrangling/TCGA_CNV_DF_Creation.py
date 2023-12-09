"""
XENABROWSER TCGA CNV DATASET IMPORT AND TRANSFORMATION

This script transforms raw TCGA CNV data into a readable CSV format,
where Genes are represented as columns and Samples as rows. It relies 
on the mRNA datasets to retrieve the same Genes of Interest (GoIs) and 
Random Control Genes (RCGs) in order to initialise the dataframe without
having to query the STRING DB.


The script follows these main steps:
1. Obtain pre-defined Set of Interest (SoI) and Control set
	from the RNA datasets that were created
2. Read and reformat raw CNV TCGA data.
3. Create the dataset focusing only on genes of interest.
4. Create the control dataset with the control set, and some additional
	random control genes in case the SoI and Control set differ in 
	gene number due to differences in the CNV and RNA datasets.
5. Convert data to numeric format and save the resulting datasets 
as CSV files.

Note: 
- The script assumes the access to the raw CNV TCGA file as well as
 the two additional Set-of-interest, and control set RNA datasets
- Output data is saved in CSV format with specific naming conventions.
"""

import random
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from io import StringIO


pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

# Loading the datasets with Pandas. 
cntrl_RNA_df = pd.read_csv('../../Data/RNA/Control_Genes_RNA_df.csv')
GoI_RNA_df = pd.read_csv('../../Data/RNA/Genes_of_Interest_RNA_Data_df.csv')


# lists containing the same genes of interest and the same
# random control genes as in the RNA dataframes.

GoI_list = list(GoI_RNA_df.columns[1:])
cntrl_list = list(cntrl_RNA_df.columns[1:])


# Read the raw CNV dataset.

with open('../../Data/CNV/TCGA_PanCan_CNV.txt', 'r') as file:
	raw_data = file.readlines()


# Process raw data to separate sample headers 
# and reformat the data into a usable structure.

samples = raw_data[0].strip('\n').split('\t')
data = raw_data[1:]
reformated_data = [item.strip('\n').split('\t') for item in data]


# List for sample CNV data of the gene-set of interest

SoI = []

for item in reformated_data:
    if item[0] in GoI_list:
        SoI.append(item)


# Construct the dataframe based on the genes of interest, 
# with genes as rows and samples as columns.
# Values are the CNVs

SoI_CNV_df = (
    pd.DataFrame(SoI, columns=samples)
	.set_index('Sample')
	.T
	.reset_index()
	.rename_axis(None, axis=1)
	.rename(columns = {'index':'Sample'})
	.drop_duplicates('Sample')
)


# List for sample CNV data from the random control genes

ctrl_set = []
leftover_ctrl_set = []


for item in reformated_data:
    if item[0] in cntrl_list:
        ctrl_set.append(item)
    elif item[0] not in GoI_list:
    	leftover_ctrl_set.append(item)


# If the number of genes in the set of interest exceeds that of
# the control-set then add a the difference from the leftover control
# list using random selection.

if len(SoI) > len(ctrl_set):
	gene_dif = len(SoI) - len(ctrl_set)
	random_cntrls = random.sample(leftover_ctrl_set, k=gene_dif)
	ctrl_set = ctrl_set + random_cntrls


# Construct the dataframe based on the random control genes, with genes 
# as rows and samples as columns. Values are the CNVs.

ctrl_set_CNV_df = (
    pd.DataFrame(ctrl_set, columns=samples)
    .drop_duplicates()
    .set_index('Sample')
    .T
    .reset_index()
    .rename_axis(None, axis=1)
    .rename(columns = {'index':'Sample'})
    .drop_duplicates('Sample')
)


# Convert all values in both dataframes to float, handling non-numeric values as NaN.

for df in [SoI_CNV_df, ctrl_set_CNV_df]:
	for col in df.columns[1:]:
		df[col] = pd.to_numeric(df[col], errors='coerce')

for col in ctrl_set_CNV_df.columns[1:]:
	ctrl_set_CNV_df[col] = pd.to_numeric(ctrl_set_CNV_df[col], errors='coerce')


print(SoI_CNV_df.iloc[:5,:5])
print(ctrl_set_CNV_df.iloc[:5,:5])


# Save the formatted gene set and control set dataframes as CSV files.

# SoI_CNV_df.to_csv('../../Data/CNV/Genes_of_Interest_CNV_df.csv', index = False)
# ctrl_set_CNV_df.to_csv('../../Data/CNV/Control_Genes_CNV_df.csv', index = False)




