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
import dask.dataframe as dd
import subprocess
import time

from STRING_SCRIPTS.STRING_GENE_LIST_function2 import Gene_list

pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)


nrows = None
interaction_score = 400

protein_genes_of_interest = [
    "MAN1","TMPO", "EMD", "TOR1AIP1", 
    'IFFO1', "LMNA", "LMNB1", "LMNB2", 
    "LBR", "LEMD2", "PRR14", "CBX5",
    "CBX1","CBX3", "LEMD3", 'BANF1',
    "BAF", "TERB1", "TERB2", "MAJN" 
    ]




# # Query STRING database to obtain a set (Set-of-Interest, SoI) of genes 
# # that interact with the products of the genes in protein_genes_of_interest.
# STRING_SoI = Gene_list(protein_genes_of_interest, interaction_score)


###############      TCGA SPECIFIC DF CREATION.       #################


# # Reading the  RNA TPM metadata dataset that contains the gene names.
# TCGA_TPM_annotations_df = pd.read_csv("../../../Data/Other/TCGA_meta/TCGA_PanCan_TPM_Annotations.csv")
# print(TCGA_TPM_annotations_df.head(),'\n')


# # Get the ENSEMBL ids that represent the Genes of interest and are in the TCGA dataset.
# SOI_ENS_ids = TCGA_TPM_annotations_df[TCGA_TPM_annotations_df['gene'].isin(STRING_SoI)].id
# print(str(SOI_ENS_ids.count()) + " TCGA genes of the total " + str(len(STRING_SoI)) + 
# 	" interacting STRING genes were found in the dataset \n")


##### TCGA full tpm dataframe


# # If the full tpm df is already created and saved, read it in.
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


## If the TCGA Set of Interest tpm df is already created and saved, read it in.
TCGA_SOI_tpm_df = pd.read_csv('../../../Data/RNA_Data/TCGA_TPM/TCGA_mRNA_TPM_SOI.csv')
print("TCGA Set of Interest tpm df",'\n',
	TCGA_SOI_tpm_df.iloc[:5,:5],'\n',
)

# # Else create the SOI tpm df.
# TCGA_SOI_tpm_df = Full_tpm_df[Full_tpm_df['Gene'].isin(STRING_SoI)].reset_index(drop = True)

# # Save the SOI Tpm df.
# TCGA_SOI_tpm_df.to_csv('../../../Data/RNA_Data/TCGA_TPM/TCGA_mRNA_TPM_SOI.csv', index = False)


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
# print(GTEx_tpm_df.shape,'\n')

# # Else create the full tpm df.
# raw_GTEx_tpm_df = pd.read_csv('../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_TPM/Raw_GTEx_RNA_TPM.txt', sep='\t', skiprows=2)
# raw_GTEx_tpm_df = raw_GTEx_tpm_df.rename(columns = {"Name":"id"})

# # Drop samples based on RIN number < 6 & SMAFRZE = "EXCLUDE"
# # Get IDs of samples to keep
# Annots_GTEX_IDs = GTEx_Annotations_df[(GTEx_Annotations_df.SMRIN > 5.9) & (GTEx_Annotations_df.SMAFRZE != "EXCLUDE")].SAMPID.to_list()

# # Drop samples based on Death hardy score = 4
# # Get IDs of samples to keep
# Meta_GTEX_SUBJIDs = GTEx_meta_df[GTEx_meta_df.DTHHRDY != 4].SUBJID.to_list()

# #Full list of IDs to keep
# GTEX_IDs = ['id','Description']

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
# GTEx_tpm_df = raw_GTEx_tpm_df.drop_duplicates('id').drop(columns = cols_to_drop).rename(columns = {"Description": "Gene"})

# # Convert the RNA values to numeric, replace 0 with 0.001, and apply log2
# GTEx_tpm_df.iloc[:, 2:] = np.log2(raw_GTEx_tpm_df.iloc[:, 2:].replace(0, 0.001).apply(pd.to_numeric, errors='coerce'))
# GTEx_tpm_df[GTEx_tpm_df.columns[2:]]= GTEx_tpm_df[GTEx_tpm_df.columns[2:]].map(lambda x: round(x, 4))

# # # Save the full tpm df.
# GTEx_tpm_df.to_csv('../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_TPM/GTEx_RNA_TPM.csv', index = False)


###### GTEx Set of interest DF


# # If the GTEx Set of Interest tpm df is already created and saved, read it in.
# GTEx_SOI_tpm_df = pd.read_csv('../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_TPM/GTEx_mRNA_TPM_SOI.csv', nrows = nrows)
# GTEx_SOI_tpm_df = GTEx_SOI_tpm_df.rename(columns = {"Name":"id"})
# print("GTEx Set of Interest tpm df",'\n',
# 	GTEx_SOI_tpm_df.iloc[:5,:5],'\n', 
# 	GTEx_SOI_tpm_df.shape,'\n'
# )

# # Else create the SOI tpm df.
# GTEx_SOI_tpm_df = GTEx_tpm_df[
# 	(GTEx_tpm_df['id'].isin(TCGA_SOI_tpm_df.id.to_list())) | 
# 	(GTEx_tpm_df['Gene'].isin(TCGA_SOI_tpm_df.Gene.to_list()))
# 	].reset_index(drop = True)

# print(GTEx_SOI_tpm_df.iloc[:10,:10])
# print(GTEx_SOI_tpm_df.shape)

# # Save the GTEx SOI Tpm df.
# GTEx_SOI_tpm_df.to_csv('../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_TPM/GTEx_mRNA_TPM_SOI.csv', index = False)


###### GTEx Full Control DF


# If the GTEx full Control tpm df is already created and saved, read it in.
GTEx_Full_Ctrl_tpm_df = pd.read_csv('../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_TPM/GTEx_mRNA_TPM_CTRL.csv', nrows = nrows)
GTEx_Full_Ctrl_tpm_df = GTEx_Full_Ctrl_tpm_df.rename(columns = {"Name":"id"})
print("GTEx Full Control tpm df",'\n',
	GTEx_Full_Ctrl_tpm_df.iloc[:5,:5],'\n', 
	GTEx_Full_Ctrl_tpm_df.shape,'\n'
)

# # Else Create the df with all control genes
# GTEx_Full_Ctrl_tpm_df = GTEx_tpm_df[~GTEx_tpm_df['Gene'].isin(GTEx_SOI_tpm_df.Gene.to_list())].reset_index(drop = True)

# # Save the Ctrl Tpm df.
# GTEx_Full_Ctrl_tpm_df.to_csv('../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_TPM/GTEx_mRNA_TPM_CTRL.csv', index = False)


###############      GTEx & TCGA SHARED CONTROLS DF CREATION.       #################


# # # If the TCGA Shared Control tpm df is already created and saved, read it in.
TCGA_Shrd_Ctrl_tpm_df = pd.read_csv('../../../Data/RNA_Data/TCGA_TPM/TCGA_mRNA_TPM_CTRL_SHRD.csv', nrows = nrows)
print("TCGA Shared Control tpm df",'\n',
	TCGA_Shrd_Ctrl_tpm_df.iloc[:5,:5],'\n',
	TCGA_Shrd_Ctrl_tpm_df.shape,'\n'
)

# # If the GTEx Shared Control tpm df is already created and saved, read it in.
GTEx_Shrd_Ctrl_tpm_df = pd.read_csv('../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_TPM/GTEx_mRNA_TPM_CTRL_SHRD.csv', nrows = nrows)
print("GTEx Shared Control tpm df",'\n',
	GTEx_Shrd_Ctrl_tpm_df.iloc[:5,:5],'\n', 
	GTEx_Shrd_Ctrl_tpm_df.shape,'\n'
)

# ## Else Create the Dataframes:
# # Identify shared genes and IDs
# shared_genes = set(TCGA_Full_Ctrl_tpm_df['Gene']).intersection(set(GTEx_Full_Ctrl_tpm_df['Gene']))
# print("Shared Genes: ", len(shared_genes), '\n')

# shared_ids = set(TCGA_Full_Ctrl_tpm_df['id']).intersection(set(GTEx_Full_Ctrl_tpm_df['id']))
# print("Shared IDs: ", len(shared_ids), '\n')

# # Filter both dataframes to only include rows with shared genes or IDs, remove any duplicate IDs 
# TCGA_Shrd_Ctrl_tpm_df = (TCGA_Full_Ctrl_tpm_df[TCGA_Full_Ctrl_tpm_df['Gene'].isin(shared_genes) | 
# 	TCGA_Full_Ctrl_tpm_df['id'].isin(shared_ids)]
# 	.drop_duplicates(subset = "id")
# 	.reset_index(drop = True)
# )

# GTEx_Shrd_Ctrl_tpm_df = (GTEx_Full_Ctrl_tpm_df[GTEx_Full_Ctrl_tpm_df['Gene'].isin(shared_genes) | 
# 	GTEx_Full_Ctrl_tpm_df['id'].isin(shared_ids)]
# 	.drop_duplicates(subset = "id")
# 	.reset_index(drop = True)
# )



# ## Save the DFs
# GTEx_Shrd_Ctrl_tpm_df.to_csv('../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_TPM/GTEx_mRNA_TPM_CTRL_SHRD.csv', index = False)
# TCGA_Shrd_Ctrl_tpm_df.to_csv('../../../Data/RNA_Data/TCGA_TPM/TCGA_mRNA_TPM_CTRL_SHRD.csv', index = False)





'''
Smaller control sets creation.
'''

# Basing the control sets lengths off of the df for the set of interest
SOI_len = len(TCGA_SOI_tpm_df.Gene.to_list())
# SOI_len = 15

all_ctrl_genes = len(TCGA_Shrd_Ctrl_tpm_df)
print("the number of all the control genes is: ", all_ctrl_genes, '\n')

# Variable to keep track of the number of control genes 
# which remain unaccounted for in any of the control datasets.
dif_ctrl = all_ctrl_genes

# Empty list to contain the gene IDs that have been used
used_Gene_ids = set()

step = 1
while dif_ctrl > 0:
	start_time = time.time()
	print('\n\n\n\n\t',"DATASET #: ", step, '\n')


	# Create a random TCGA control df the same size as the SOI using all the shared control genes 
	rndm_TCGA_ctrl_set = TCGA_Shrd_Ctrl_tpm_df.sample(n = SOI_len).reset_index(drop = True)
	len_rndm_TCGA_ctrl_set = len(rndm_TCGA_ctrl_set)
	print("Random Tcga Control Set",'\n', 
		# rndm_TCGA_ctrl_set.iloc[:5,:5], '\n',
		"Gene count: ", len_rndm_TCGA_ctrl_set, '\n'
	)

	# Variable to keep track of the genes that are being assigned to the new df
	rndm_ids = rndm_TCGA_ctrl_set.id.to_list()

	# Variable to keep track of the number of newly added ids
	new_ids = len(set([ID for ID in rndm_ids if ID not in used_Gene_ids]))

	# If all the genes in the new random df have already been used, 
	# then start over.
	if new_ids == 0:
		print("No new genes were found", '\n')

		end_time = time.time()
		duration = end_time - start_time
		print("Duration:", duration, "seconds",'\n\n\n')
		continue

	# If some genes have not been used, 
	# then proceed
	else:
		print("The number of new genes in this set is: ", new_ids,'\n')

		# Create the corresponding GTEx control df based on the IDs
		Corres_GTEx_ctrl_set = GTEx_Shrd_Ctrl_tpm_df[GTEx_Shrd_Ctrl_tpm_df.id.isin(rndm_ids)].reset_index(drop = True)
		len_Corres_GTEx_ctrl_set = len(Corres_GTEx_ctrl_set)
		print("Corresponding GTEx control set",'\n',
			# Corres_GTEx_ctrl_set.iloc[:5,:5], '\n',
			len_Corres_GTEx_ctrl_set, '\n'
		)

		# If there is a mismatch between the TCGA and GTEx controls,
		# fill in the difference using the missing gene names instead of IDs
		if len_rndm_TCGA_ctrl_set != len_Corres_GTEx_ctrl_set:
			print("Mistmatch identified", '\n')
			remaining_genes = list(rndm_TCGA_ctrl_set[~rndm_TCGA_ctrl_set.id.isin(Corres_GTEx_ctrl_set.id.to_list())].Gene)

			# Number of genes that need to be added to the control df
			remaining_number = len(remaining_genes)
			print(" Remaining genes to be filled :", '\n',
			 # remaining_genes,'\n' ,
			 "Number: ",remaining_number, '\n'
			)

			# Dataframe with the data from the remaining identified genes
			remaining_genes_df = GTEx_Shrd_Ctrl_tpm_df[GTEx_Shrd_Ctrl_tpm_df.Gene.isin(remaining_genes)]
			# print(" Remaining genes dataframe", '\n',
			#  remaining_genes_df.iloc[:5,:5], '\n',
			#  len(remaining_genes_df), '\n'
			# )

			# If the number of genes in the dataframe built using the remaining genes 
			# is smaller than the required number of genes, then start over.
			if len(remaining_genes_df) < remaining_number:
				print("Not enough genes in the remaining genes df to fill in the control", '\n')

				end_time = time.time()
				duration = end_time - start_time
				print("Duration:", duration, "seconds",'\n\n\n')
				continue

			else:
				# Subset of the remaining df with the correct number of genes
				subset_remaining_genes_df = remaining_genes_df.sample(n = remaining_number)
				# print(" Subset of Remaining genes dataframe", '\n',
				#  subset_remaining_genes_df.iloc[:5,:5], '\n',
				#  len(subset_remaining_genes_df), '\n'
				# )

				Filled_Corres_GTEx_ctrl_set = pd.concat([Corres_GTEx_ctrl_set, subset_remaining_genes_df]).reset_index(drop = True)
				print(" Filled dataframe", '\n',
				 # Filled_Corres_GTEx_ctrl_set.iloc[:5,:5], '\n',
				 "Gene count: ",len(Filled_Corres_GTEx_ctrl_set), '\n'
				)

				# Add the IDs to the list of used gene ids
				All_used = Filled_Corres_GTEx_ctrl_set.id.to_list() + rndm_TCGA_ctrl_set.id.to_list()
				used_Gene_ids.update(All_used)

				#Save the dfs
				rndm_TCGA_ctrl_set.to_csv(f'../../../Data/RNA_Data/TCGA_TPM/TCGA_control_sets/TCGA_TPM_RNA_Control_df{step}.csv', index = False)
				Filled_Corres_GTEx_ctrl_set.to_csv(f'../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_TPM/GTEx_control_sets/GTEx_TPM_RNA_Control_df{step}.csv', index = False)




		# If the TCGA and GTEx controls matchup, save the sets
		else: 
			print("All Good",'\n')
			rndm_TCGA_ctrl_set.to_csv(f'../../../Data/RNA_Data/TCGA_TPM/TCGA_control_sets/TCGA_TPM_RNA_Control_df{step}.csv', index = False)
			Corres_GTEx_ctrl_set.to_csv(f'../../../Data/RNA_Data/GTEx_RNA/GTEx_RNA_TPM/GTEx_control_sets/GTEx_TPM_RNA_Control_df{step}.csv', index = False)

		print("Number of genes that have been used: ", len(used_Gene_ids),'\n')

		# Update the difference variable
		dif_ctrl = all_ctrl_genes - len(used_Gene_ids)
		print("The Number of genes that have yet to be used: ", dif_ctrl,'\n')

		end_time = time.time()
		duration = end_time - start_time
		print("Duration:", duration, "seconds",'\n\n\n')
		step +=1 




# return_code = subprocess.call(["afplay", "/Users/Dyll/Downloads/mixkit-hawk-bird-squawk-1268.wav"])


