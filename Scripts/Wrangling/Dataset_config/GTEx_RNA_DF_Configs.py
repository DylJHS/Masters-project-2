import random
import pandas as pd
import pyarrow as pa
from io import StringIO


pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

"""
GTEx Configuration Script

The purpose of the script is to clean the GTEx RNA dataset so that
it contains only samples from heathly subjects defined by the GTEx Portal,
and to create and save the respective datasets that can be used as healthy controls. 
This includes:
- A dataset based solely on the Genes of Interest
- Multiple controls whose number is the same as the Cancerous ones
and whose composition is as close to the TCGA controls as possible. 

"""

# Reading the full RNA GTEx dataset
GTEx_Samp_df = pd.read_csv("../../../Data/RNA_Data/GTEx_RNA/Raw_GTEx_RNA_Norm.csv", index_col = 0)

# Reading the first RNA GTEx metadata dataset that contains the phenotype of the subject
GTEx_meta_df = pd.read_csv("../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_Subject_meta.csv")

# Reading the second metadata dataset that contains the rest of the metadata for the subjects
GTEx_Annotations_df = pd.read_csv("../../../Data/RNA_Data/GTEx_RNA/GTEx_Annotations.csv")

### The GTEx dataset contains 17384 columns (subjects) ###




"""
Gathering the information that is needed from the various datasets 
"""

# Get the SUBJIDS that correspond to the samples to be used
Healthy_Subjects = GTEx_meta_df[GTEx_meta_df['DTHHRDY'].isin([1, 2, 3, 0])]['SUBJID'].unique()

# Get the column names from the full RNA GTEX dataset
cols = GTEx_Samp_df.columns[3:]

# Get the columns that should be dropped from the dataset
cols_to_drop = [col for col in cols if not col[:10] in Healthy_Subjects]

# Create a new DataFrame without the columns to be dropped
GTEx_Healthy_df = GTEx_Samp_df.drop(columns=cols_to_drop)

### The GTEx Healthy dataframe contains 10696 columns (subjects) ###



# Get the Sample IDs that correspond to the samples to be used
selected_samples = GTEx_Annotations_df[GTEx_Annotations_df['SMRIN'] >= 6]['SAMPID']
selected_samples_list = ["Description"] + selected_samples.tolist()

# Create a DataFrame with the columns to be used
GTEx_Viable_Healthy = GTEx_Healthy_df[[col for col in selected_samples_list if col in GTEx_Healthy_df]]
### The GTEx Healthy dataframe contains 9645 columns (subjects) ###


GTEx_Viable_Healthy = GTEx_Viable_Healthy.groupby('Description').mean().reset_index().set_index('Description')
GTEx_Viable_Healthy.index = GTEx_Viable_Healthy.index.str.upper()

# print(GTEx_Viable_Healthy.iloc[:5,:5])
# print(GTEx_Viable_Healthy.shape)



# Get the genes that are used in the TCGA PanCan Dataset

TCGA_Pancan_df =  pd.read_csv("../../../Data/RNA_Data/TCGA_mRNA_Norm.csv", usecols = [0])
### 17885 Genes are common to both datasets, with the TCGA set containing a total of 20530 unique genes ###


# Get the Genes of Interest that are used in the TCGA PanCan 

SoI_list =  pd.read_csv('../../../Data/RNA_Data/Genes_of_Interest_RNA_Data_df.csv', index_col=0, nrows=0).columns.tolist()
SoI_list = [word.upper() for word in SoI_list]
### 304 Genes are common to both the TCGA SoI and the GTEx gene-set ###


# Create the dataframe for the GTEx using the genes of interest.

GTEx_SOI_df = (GTEx_Viable_Healthy[GTEx_Viable_Healthy.index.isin(SoI_list)]
	.T
	.reset_index(drop = False)
	.rename(columns = {'index':'SampleID'})
	.rename_axis(None, axis = 1)
	)

for col in GTEx_SOI_df.columns[1:]:
	GTEx_SOI_df[col] = pd.to_numeric(GTEx_SOI_df[col], errors='coerce')

# GTEx_SOI_df.to_csv('../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_Healthy/GTEx_RNA_SoI.csv', index = False)

# print(GTEx_SOI_df.iloc[:5,:5])
# print(GTEx_SOI_df.shape)



# Create the dataframe for the GTEx using the genes that are not of interest.

GTEx_Cntrl_df = (GTEx_Viable_Healthy[~GTEx_Viable_Healthy.index.isin(SoI_list)]
	.T
	.reset_index(drop = False)
	.rename(columns = {'index':'SampleID'})
	.rename_axis(None, axis = 1)
	)

for col in GTEx_Cntrl_df.columns[1:]:
	GTEx_Cntrl_df[col] = pd.to_numeric(GTEx_Cntrl_df[col], errors='coerce')

# GTEx_Cntrl_df.to_csv('../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_Healthy/Controls/GTEx_RNA_Full_Cntrl.csv', index = False)

print(GTEx_SOI_df.iloc[:5,:5])
print(GTEx_SOI_df.shape)


# # Create the control sets that match up as closely as possible to the TCGA controls

# If the Full control set has already been created
GTEx_Cntrl_df = pd.read_csv('../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_Healthy/Controls/GTEx_RNA_Full_Cntrl.csv')

n = 1

while n < 4:

	try:

		# Load the mRNA control dataset. 
		# list containing the same genes as the RNA control df.

		TCGA_cntrl =  pd.read_csv(f'../../../Data/RNA_Data/Control_Gene_Sets/Control_Genes_RNA_Data_df{n}.csv', index_col=0, nrows=0).columns.tolist()
		TCGA_cntrl_list =  ['SampleID'] + TCGA_cntrl



		GTEx_ctrl_set = GTEx_Cntrl_df[[col for col in TCGA_cntrl_list if col in GTEx_Cntrl_df]]
		GTEx_ctrl_set_list = list(GTEx_ctrl_set.columns)


		leftover_ctrl_set = GTEx_Cntrl_df[[col for col in GTEx_Cntrl_df if col not in TCGA_cntrl_list]]
		leftover_ctrl_set_list = list(leftover_ctrl_set.columns)


		# If the number of genes in the TCGA control exceeds that of
		# the GTEx control-set then add the difference from the leftover control
		# list using random selection.

		gene_dif = len(TCGA_cntrl_list) - len(GTEx_ctrl_set_list)

		print(" The number of genes that are commmon to both TCGA and GTEx sets: ", gene_dif)

		if gene_dif > 0:

			random_genes = random.sample(leftover_ctrl_set_list, k=gene_dif)
			all_ctrls = list(GTEx_ctrl_set_list + random_genes)

			GTEx_ctrl_new_set = GTEx_Cntrl_df[all_ctrls]

			print(GTEx_ctrl_new_set.iloc[:5,:5], '\n')


			# Save the formatted gene set and control set dataframes as CSV files.

			# GTEx_ctrl_new_set.to_csv(f'../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_Healthy/Controls/GTEx_Control_RNA_Data_df{n}.csv', index = False)

		n += 1

	except FileNotFoundError:
		print("\n\n\n Last dataset created \n\n\n")
		break



