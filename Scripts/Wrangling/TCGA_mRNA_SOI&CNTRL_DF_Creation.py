"""
XENABROWSER TCGA DATASET IMPORT AND TRANSFORMATION

This script transforms raw TCGA data (RNAseq, CNV) into a 
Set-of-Interest mRNA dataset, and multiple control mRNA datasets. 
The datasets are maintained in a CSV format where genes are represented 
as columns and samples as rows. 
It leverages the Gene_list function to retrieve a list of genes whose 
products interact with a specified set of proteins, using the STRING database 
for querying.

The script follows these main steps:
1. Query STRING database for interacting genes based on 
	a predefined list of protein genes.
2. Read and reformat raw TCGA data.
3. Create a dataset focusing only on genes of interest.
4. Generate control datasets with a random selection of genes, ensuring that each
	matches the shape of the gene of interest dataset (SoI_df). 
	Control datasets are generated until each of the genes in the complete
	TCGA file that are not part of the SoI are included at least once in a 
	control dataset.
5. Convert data to numeric format and save the resulting datasets as CSV files.

Note: 
- The script assumes the availability of raw TCGA data and requires an internet connection
  for querying the STRING database.
- Output data is saved in CSV format with specific naming conventions.
"""

import random
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from io import StringIO

from STRING_GENE_LIST import Gene_list

pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

Data_type = 'RNA_Data'



protein_genes_of_interest = [
    "MAN1","TMPO", "EMD", "TOR1AIP1", 
    'IFFO1', "LMNA", "LMNB1", "LMNB2", 
    "LBR", "LEMD2", "PRR14", "CBX5"
]


# Query STRING database to obtain a set (Set-of-Interest, SoI) of genes 
# that interact with the products of the genes in protein_genes_of_interest.

STRING_SoI = Gene_list(protein_genes_of_interest)


# Read raw TCGA mRNA data from a file and store each line in raw_data.

with open('../../Data/RNA_Data/TCGA_mRNA_Norm', 'r') as file:
	raw_data = file.readlines()


# Process raw data to separate sample headers 
# and reformat the data into a usable structure.

samples = raw_data[0].strip('\n').split('\t')
data = raw_data[1:]
reformated_data = [item.strip('\n').split('\t') for item in data]
# print(len(reformated_data))


# Create a list of data for genes of interest, 
# ensuring each gene is unique.

SoI_seen = set()
SoI = []

for item in reformated_data:
    if item[0] in STRING_SoI and item[0] not in SoI_seen:
        SoI.append(item)
        SoI_seen.add(item[0])


# Construct a dataframe from genes of interest, 
# with genes as rows and samples as columns.

SoI_df = (
    pd.DataFrame(SoI, columns=samples)
	.set_index('sample')
	.T
	.reset_index()
	.rename_axis(None, axis=1)
	.rename(columns = {'index':'Sample'})
	.drop_duplicates('Sample')
)


# Convert all values in dataframe to float, handling non-numeric values as NaN.

for col in SoI_df.columns[1:]:
	SoI_df[col] = pd.to_numeric(SoI_df[col], errors='coerce')


# Save the formatted gene set and control set dataframes as CSV files.

# SoI_df.to_csv(f'../../Data/{Data_type}/Genes_of_Interest_{Data_type}_df.csv', index = False)



# Generate a control list containing the data for the 
# genes not in genes of interest.

ctrl_set_seen = set()
ctrl_set = []

for gene_data in reformated_data:
    if gene_data[0] not in STRING_SoI and gene_data[0] not in ctrl_set_seen:
        ctrl_set.append(gene_data)
        ctrl_set_seen.add(gene_data[0])


# print(len(ctrl_set_seen), len(ctrl_set))



# Create a list containing the control genes that have
# already appeared in a control dataset.

ctrl_used = set()


# Variable to keep track of the number of control genes 
# which remain unaccounted for in any of the control datasets.

dif_ctrl = len(ctrl_set_seen) - len(ctrl_used)


# Variable to name and keep track of the current control dataset

n=1


# Control datasets are created while there are still control genes 
# that have not yet been used in a control dataset.

while dif_ctrl != 0:

	# Randomly select genes for the control set to match the 
	# number of genes in the SoI. Repeat selection until the control
	# set dataframe matches the shape of the gene set dataframe.

	while True:

		rndm_ctrl_set = random.sample(ctrl_set, k=len(SoI))


		# Variable to keep track of the number of genes in this control set 
		# which do not show up in any other control dataset

		num_discovered_genes = len(
			[gene[0] for gene in rndm_ctrl_set if gene[0] not in ctrl_used]
		)


		# Do not create control datasets that do not contain a single 
		# previously unseen gene
		
		if num_discovered_genes == 0:
			continue

		var = f'rndm_ctrl_set_df{n}'

		globals()[var] = (
		    pd.DataFrame(rndm_ctrl_set, columns=samples)
		    .drop_duplicates()
		    .set_index('sample')
		    .T
		    .reset_index()
		    .rename_axis(None, axis=1)
		    .rename(columns = {'index':'Sample'})
		    .drop_duplicates('Sample')
		)


		# Convert all values in both dataframes to float, handling non-numeric values as NaN.

		for col in globals()[var].columns[1:]:
			globals()[var][col] = pd.to_numeric(globals()[var][col], errors='coerce')
			# print(globals()[var].shape, globals()[var].Sample.nunique()) 



		if globals()[var].shape == SoI_df.shape:
			print('\n'
				dif_ctrl, 
				'\n',
				num_discovered_genes, 
				'\n', 
				n, 
				'\n')

			print(globals()[var].iloc[:3,:10])

			for gene_data in rndm_ctrl_set:
						ctrl_used.add(gene_data[0])

			dif_ctrl = len(ctrl_set_seen) - len(ctrl_used)

			# Save the formatted gene set and control set dataframes as CSV files.

			# globals()[var].to_csv(f'../../Data/{Data_type}/Control_Genes_{Data_type}_df{n}.csv', index = False)
			n += 1
			break  # Exit the loop if the shapes are equal



