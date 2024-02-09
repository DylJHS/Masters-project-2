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
import numpy as np

from STRING_SCRIPTS.STRING_GENE_LIST_function2 import Gene_list

pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

interaction_score = 400

protein_genes_of_interest = [
    "MAN1","TMPO", "EMD", "TOR1AIP1", 
    'IFFO1', "LMNA", "LMNB1", "LMNB2", 
    "LBR", "LEMD2", "PRR14", "CBX5",
    "CBX1","CBX3", "LEMD3"
    ]


nrows = None

# Query STRING database to obtain a set (Set-of-Interest, SoI) of genes 
# that interact with the products of the genes in protein_genes_of_interest.
STRING_SoI = Gene_list(protein_genes_of_interest, interaction_score)


###############      TCGA SPECIFIC DF CREATION.       #################


# Reading the  RNA TPM metadata dataset that contains the gene names.
TCGA_TPM_annotations_df = pd.read_csv("../../../Data/Other/TCGA_meta/TCGA_PanCan_TPM_Annotations.csv")
# print(TCGA_TPM_annotations_df.head(),'\n')


# # Get the ENSEMBL ids that represent the Genes of interest and are in the TCGA dataset.
# SOI_ENS_ids = TCGA_TPM_annotations_df[TCGA_TPM_annotations_df['gene'].isin(STRING_SoI)].id
# print(str(SOI_ENS_ids.count()) + " TCGA genes of the total " + str(len(STRING_SoI)) + 
# 	" interacting STRING genes were found in the dataset \n")


##### TCGA full tpm dataframe


## If the full tpm df is already created and saved, read it in.
# Full_tpm_df = pd.read_csv('../../../Data/RNA_Data/TCGA_TPM/TCGA_mRNA_TPM.csv', nrows = nrows)
# print(Full_tpm_df.iloc[:5,:5],'\n')
# print(Full_tpm_df.shape,'\n')

# # ## Else create the full tpm df.
# raw_tpm_df = pd.read_csv('../../../Data/RNA_Data/TCGA_TPM/tcga_RSEM_gene_tpm.txt', sep = "\t")

# # list of the tumour types coodes that are NOT allowed in the analysis 
# tumour_type_codes = ['11','15','16','20','50','60','61','99']

# # list of the TSS codes that are NOT allowed in the analysis 
# tss_codes = ['07','AV']

# # Check if the last 6th & 7th digits are in the list to drop
# tss_codes_to_drop = raw_tpm_df.columns[raw_tpm_df.columns.str[5:7].isin(tss_codes)].to_list()

# # Check if the last two digits are  in the list to drop
# tumour_type_codes_to_drop = raw_tpm_df.columns[raw_tpm_df.columns.str[-2:].isin(tumour_type_codes)].to_list()

# #combine the columns to drop
# columns_to_drop = tss_codes_to_drop + tumour_type_codes_to_drop
# columns_to_drop.append('sample')

# Full_tpm_df = (TCGA_TPM_annotations_df[['id','gene']]
# 	.merge(raw_tpm_df, left_on = "id", right_on = "sample")
#   .drop_duplicates('sample')
# 	.drop(columns = columns_to_drop)
# 	.rename(columns = {"gene": "Gene"})
# 	.reset_index(drop = True)
# 	)
# # print(Full_tpm_df.iloc[:5,:5],'\n')
# # print(Full_tpm_df.shape,'\n')

# # ## Save the Full Tpm df.
# Full_tpm_df.to_csv('../../../Data/RNA_Data/TCGA_TPM/TCGA_mRNA_TPM.csv', index = False)


##### TCGA Set of Interest tpm dataframe


# ## If the TCGA Set of Interest tpm df is already created and saved, read it in.
# SOI_tpm_df = pd.read_csv('../../../Data/RNA_Data/TCGA_TPM/TCGA_mRNA_TPM_SOI.csv')
# print("TCGA Set of Interest tpm df",'\n',SOI_tpm_df.iloc[:5,:5],'\n')
# print(SOI_tpm_df.shape,'\n')

# # Else create the SOI tpm df.
# SOI_tpm_df = Full_tpm_df[Full_tpm_df['Gene'].isin(STRING_SoI)].reset_index(drop = True)

# # Save the SOI Tpm df.
# SOI_tpm_df.to_csv('../../../Data/RNA_Data/TCGA_TPM/TCGA_mRNA_TPM_SOI.csv', index = False)


##### TCGA full Control tpm dataframe


# # If the TCGA full Control tpm df is already created and saved, read it in.
TCGA_Full_Ctrl_tpm_df = pd.read_csv('../../../Data/RNA_Data/TCGA_TPM/TCGA_mRNA_TPM_CTRL.csv', nrows = nrows)
print("TCGA full Control tpm df",'\n', TCGA_Full_Ctrl_tpm_df.iloc[:5,:5],'\n')
print(TCGA_Full_Ctrl_tpm_df.shape,'\n')

# # Else Create the df with all control genes
# TCGA_Full_Ctrl_tpm_df = Full_tpm_df[~Full_tpm_df['Gene'].isin(STRING_SoI)].reset_index(drop = True)

# # Save the Ctrl Tpm df.
# TCGA_Full_Ctrl_tpm_df.to_csv('../../../Data/RNA_Data/TCGA_TPM/TCGA_mRNA_TPM_CTRL.csv', index = False)


###############      GTEx SPECIFIC DF CREATION.       #################


# # Reading the  RNA GTEx metadata dataset that contains the gene names.
# GTEx_meta_df = pd.read_csv("../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_Subject_meta.csv")
# GTEx_Annotations_df = pd.read_csv("../../../Data/RNA_Data/GTEx_RNA/GTEx_Annotations.csv")


###### GTEx Full TPM DF 


# ## If the full GTEx tpm df is already created and saved, read it in.
# GTEx_tpm_df = pd.read_csv('../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_TPM/GTEx_RNA_TPM.csv', nrows = nrows)
# print("GTEx full tpm df",'\n',GTEx_tpm_df.iloc[:5,:5],'\n')
# # print(GTEx_tpm_df.shape,'\n')

# # Else create the full tpm df.
# raw_GTEx_tpm_df = pd.read_csv('../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_TPM/Raw_GTEx_RNA_TPM.txt', sep='\t', skiprows=2)

# # Drop samples based on RIN number < 6 & SMAFRZE = "EXCLUDE"
# # Get IDs of samples to keep
# Annots_GTEX_IDs = GTEx_Annotations_df[(GTEx_Annotations_df.SMRIN > 5.9) & (GTEx_Annotations_df.SMAFRZE != "EXCLUDE")].SAMPID.to_list()

# # Drop samples based on Death hardy score = 4
# # Get IDs of samples to keep
# Meta_GTEX_SUBJIDs = GTEx_meta_df[GTEx_meta_df.DTHHRDY != 4].SUBJID.to_list()

# #Full list of IDs to keep
# GTEX_IDs = ['Name','Description']

# # Iterate over each item in the Meta_GTEX_SUBJIDs list
# for SUBJID in Meta_GTEX_SUBJIDs:
#     # Check if the code appears as a substring in any column name
#     matching_columns = [col for col in Annots_GTEX_IDs if SUBJID in col]
#     # Extend the GTEX_IDs list with the matching column names
#     GTEX_IDs.extend(matching_columns)

# # Convert GTEX_IDs to a set to remove duplicates
# GTEX_IDs = list(set(GTEX_IDs))

# # Get list of the IDs to remove
# cols_to_drop = [col for col in raw_GTEx_tpm_df.columns if col not in GTEX_IDs]

# # Drop the columns from the DF
# GTEx_tpm_df = raw_GTEx_tpm_df.drop_duplicates('Name').drop(columns = cols_to_drop).rename(columns = {"Description": "Gene"})

# # Convert the RNA values to numeric, replace 0 with 0.001, and apply log2
# GTEx_tpm_df.iloc[:, 2:] = np.log2(raw_GTEx_tpm_df.iloc[:, 2:].replace(0, 0.001).apply(pd.to_numeric, errors='coerce'))

# # Save the full tpm df.
# GTEx_tpm_df.to_csv('../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_TPM/GTEx_RNA_TPM.csv', index = False)


###### GTEx Set of interest DF


## If the GTEx Set of Interest tpm df is already created and saved, read it in.
# GTEx_SOI_tpm_df = pd.read_csv('../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_TPM/GTEx_mRNA_TPM_SOI.csv')
# print("GTEx Set of Interest tpm df",'\n',GTEx_SOI_tpm_df.iloc[:5,:5],'\n')
# print(GTEx_SOI_tpm_df.shape,'\n')

# # Else create the SOI tpm df.
# GTEx_SOI_tpm_df = GTEx_tpm_df[
# 	(GTEx_tpm_df['Name'].isin(SOI_tpm_df.id.to_list())) | 
# 	(GTEx_tpm_df['Gene'].isin(SOI_tpm_df.Gene.to_list()))
# 	].reset_index(drop = True)

# print(GTEx_SOI_tpm_df.iloc[:10,:10])
# print(GTEx_SOI_tpm_df.shape)


# # Save the GTEx SOI Tpm df.
# GTEx_SOI_tpm_df.to_csv('../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_TPM/GTEx_mRNA_TPM_SOI.csv', index = False)


###### GTEx Full Control DF


# # If the GTEx full Control tpm df is already created and saved, read it in.
GTEx_Full_Ctrl_tpm_df = pd.read_csv('../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_TPM/GTEx_mRNA_TPM_CTRL.csv', nrows = nrows)
print("GTEx Full Control tpm df",'\n',GTEx_Full_Ctrl_tpm_df.iloc[:5,:5],'\n')
print(GTEx_Full_Ctrl_tpm_df.shape,'\n')

# # Else Create the df with all control genes
# GTEx_Full_Ctrl_tpm_df = GTEx_tpm_df[~GTEx_tpm_df['Gene'].isin(GTEx_SOI_tpm_df.Gene.to_list())].reset_index(drop = True)


# # Save the Ctrl Tpm df.
# GTEx_Full_Ctrl_tpm_df.to_csv('../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_TPM/GTEx_mRNA_TPM_CTRL.csv', index = False)
# print(GTEx_Full_Ctrl_tpm_df.iloc[:5,:5],'\n')
# print(GTEx_Full_Ctrl_tpm_df.shape,'\n')




###############      GTEx & TCGA SHARED CONTROLS DF CREATION.       #################


# # # If the TCGA Shared Control tpm df is already created and saved, read it in.
# TCGA_Shrd_Ctrl_tpm_df = pd.read_csv('../../../Data/RNA_Data/TCGA_TPM/TCGA_mRNA_TPM_CTRL_SHRD.csv', nrows = nrows)
# print("TCGA Shared Control tpm df",'\n',TCGA_Shrd_Ctrl_tpm_df.iloc[:5,:5],'\n')
# print(TCGA_Shrd_Ctrl_tpm_df.shape,'\n')


# # # If the GTEx Shared Control tpm df is already created and saved, read it in.
# GTEx_Shrd_Ctrl_tpm_df = pd.read_csv('../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_TPM/GTEx_mRNA_TPM_CTRL_SHRD.csv', nrows = nrows)
# print("GTEx Shared Control tpm df",'\n',GTEx_Shrd_Ctrl_tpm_df.iloc[:5,:5],'\n')
# print(GTEx_Shrd_Ctrl_tpm_df.shape,'\n')

# Identify shared genes and IDs
shared_genes = set(TCGA_Full_Ctrl_tpm_df['Gene']).intersection(set(GTEx_Full_Ctrl_tpm_df['Gene']))
shared_ids = set(TCGA_Full_Ctrl_tpm_df['ID']).intersection(set(GTEx_Full_Ctrl_tpm_df['ID']))

# # ELSE
# # Merge the control dfs based on the Genes
# merge_gene_df = pd.merge(
# 	TCGA_Full_Ctrl_tpm_df.iloc[:,:2], 
# 	GTEx_Full_Ctrl_tpm_df.iloc[:,:2], 
# 	on='Gene', 
# 	how='inner'
# 	).drop(columns = "Name")
# print("Gene Merge DF",'\n', merge_gene_df.iloc[:5,:5],'\n')
# print(merge_gene_df.shape,'\n')

# # Merge the control dfs based on the IDs
# merge_id_df = pd.merge(
#     TCGA_Full_Ctrl_tpm_df.iloc[:, :2], 
#     GTEx_Full_Ctrl_tpm_df.iloc[:, :2], 
#     left_on='id', 
#     right_on='Name', 
#     how='inner'
# ).drop(columns=["Name", "Gene_y"]).rename(columns={"Gene_x": "Gene"})
# print("ID Merge DF",'\n', merge_id_df.iloc[:5,:5],'\n')
# print(merge_id_df.shape,'\n')

# # Concatenate the merged DataFrames and remove duplicates
# merged_df = pd.concat([merge_gene_df, merge_id_df]).drop_duplicates().reset_index(drop=True)
# print("Concatenated DF",'\n', merged_df.iloc[:5,:5],'\n')
# print(merged_df.shape,'\n')

# ## Create the Shared TCGA Control TPM df
# TCGA_Shrd_Ctrl_tpm_df = pd.merge(
# 	merged_df, TCGA_Full_Ctrl_tpm_df, 
# 	left_on = ["id", "Gene"], 
# 	right_on = ["id", "Gene"], 
# 	how='inner'
# 	).drop_duplicates(subset="Gene"
# 	).drop_duplicates(subset="id"
# )

# print("Shared TCGA Control TPM df",'\n', TCGA_Shrd_Ctrl_tpm_df.iloc[:5,:5],'\n')
# print(TCGA_Shrd_Ctrl_tpm_df.shape,'\n')


# ## Create the Shared GTEx Control TPM df
# GTEx_Shrd_Ctrl_tpm_df = pd.merge(
# 	merged_df, GTEx_Full_Ctrl_tpm_df, 
# 	left_on = ["id", "Gene"], 
# 	right_on = ["Name", "Gene"], 
# 	how='inner'
# 	).drop_duplicates(subset="Gene"
# 	).drop_duplicates(subset="id"
# 	).drop(columns = "Name"
# )


# print("Shared GTEx Control TPM df",'\n', GTEx_Shrd_Ctrl_tpm_df.iloc[:5,:5],'\n')
# print(GTEx_Shrd_Ctrl_tpm_df.shape,'\n')




# TCGA_Shrd_Ctrl_tpm_df = TCGA_Full_Ctrl_tpm_df[
# 	(TCGA_Full_Ctrl_tpm_df['Gene'].isin(GTEx_Full_Ctrl_tpm_df.Gene.to_list())) |
# 	(TCGA_Full_Ctrl_tpm_df['id'].isin(GTEx_Full_Ctrl_tpm_df.Name.to_list()))
# 	].reset_index(drop = True)

# # # # Save the TCGA Shared Ctrl Tpm df.
# # TCGA_Shrd_Ctrl_tpm_df.to_csv('../../../Data/RNA_Data/TCGA_TPM/TCGA_mRNA_TPM_CTRL_SHRD.csv', index = False)
# # print(TCGA_Shrd_Ctrl_tpm_df.iloc[:5,:5],'\n')
# # print(TCGA_Shrd_Ctrl_tpm_df.shape,'\n')


# ##### GTEx Shared Control tpm dataframe (shared genes with the TCGA control df)




# # # # Else Create the GTEx Shared Ctrl Tpm df.
# # GTEx_Shrd_Ctrl_tpm_df = GTEx_Full_Ctrl_tpm_df[
# # 	(GTEx_Full_Ctrl_tpm_df['Gene'].isin(TCGA_Full_Ctrl_tpm_df.Gene.to_list())) |
# # 	(GTEx_Full_Ctrl_tpm_df['Name'].isin(TCGA_Full_Ctrl_tpm_df.id.to_list()))
# # 	].reset_index(drop = True)

# # # # Save the GTEx Shared Ctrl Tpm df.
# # GTEx_Shrd_Ctrl_tpm_df.to_csv('../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_TPM/GTEx_mRNA_TPM_CTRL_SHRD.csv', index = False)
# # print(GTEx_Shrd_Ctrl_tpm_df.iloc[:5,:5],'\n')
# # print(GTEx_Shrd_Ctrl_tpm_df.shape,'\n')






# '''
# Smaller control sets creation.
# '''

# # # List containing the control genes that have
# # # already appeared in a control dataset.
# # ctrl_used = set()


# # # Variable to keep track of the number of control genes 
# # # which remain unaccounted for in any of the control datasets.
# # dif_ctrl = len(ctrl_set_seen) - len(ctrl_used)


# # # Variable to name and keep track of the current control dataset

# # n=1


# # # Control datasets are created while there are still control genes 
# # # that have not yet been used in a control dataset.

# # while dif_ctrl != 0:

# # 	# Randomly select genes for the control set to match the 
# # 	# number of genes in the SoI. Repeat selection until the control
# # 	# set dataframe matches the shape of the gene set dataframe.

# # 	while True:

# # 		rndm_ctrl_set = random.sample(ctrl_set, k=len(SoI))


# # 		# Variable to keep track of the number of genes in this control set 
# # 		# which do not show up in any other control dataset

# # 		num_discovered_genes = len(
# # 			[gene[0] for gene in rndm_ctrl_set if gene[0] not in ctrl_used]
# # 		)


# # 		# Do not create control datasets that do not contain a single 
# # 		# previously unseen gene
		
# # 		if num_discovered_genes == 0:
# # 			continue

# # 		var = f'rndm_ctrl_set_df{n}'

# # 		globals()[var] = (
# # 		    pd.DataFrame(rndm_ctrl_set, columns=samples)
# # 		    .drop_duplicates()
# # 		    .set_index('sample')
# # 		    .T
# # 		    .reset_index()
# # 		    .rename_axis(None, axis=1)
# # 		    .rename(columns = {'index':'Sample'})
# # 		    .drop_duplicates('Sample')
# # 		)


# # 		# Convert all values in both dataframes to float, handling non-numeric values as NaN.

# # 		for col in globals()[var].columns[1:]:
# # 			globals()[var][col] = pd.to_numeric(globals()[var][col], errors='coerce')
# # 			# print(globals()[var].shape, globals()[var].Sample.nunique()) 



# # 		if globals()[var].shape == SoI_df.shape:
# # 			print('\n',
# # 				'Iteration number: ',n,
# # 				'\n',
# # 				'number of unused genes in the full control set that remain:', dif_ctrl, 
# # 				'\n',
# # 				'number of genes in the new control set which are not part of any previous control set:',num_discovered_genes, 
# # 				'\n', 
# # 				'Dimensions of the new control set: ', SoI_df.shape,
# # 				'\n')

# # 			print(globals()[var].iloc[:3,:10])

# # 			for gene_data in rndm_ctrl_set:
# # 						ctrl_used.add(gene_data[0])

# # 			dif_ctrl = len(ctrl_set_seen) - len(ctrl_used)

# # 			print(dif_ctrl)

# # 			# Save the formatted gene set and control set dataframes as CSV files.

# # 			# globals()[var].to_csv(f'../../Data/{Data_type}/Control_Genes_{Data_type}_df{n}.csv', index = False)
# # 			n += 1
# # 			break  # Exit the loop if the shapes are equal



