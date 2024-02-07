"""
This script performs a Dataframe creation. 
Taking the raw GTEx RNA dataset and configuring it to contain 
all the genes as rows and the samples as columns.

"""

import pandas as pd
from io import StringIO


pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)


# Read raw CNV data from a file and store each line in raw_data.

with open('../../../Data/RNA_Data/tcga_RSEM_gene_tpm', 'r') as file:
	raw_data = file.readlines()



# Process raw data to separate sample headers 
# and reformat the data into a usable structure.

columns = raw_data[0].strip('\n').split('\t')
values = raw_data[1:]
reformated_values = [item.strip('\n').split('\t') for item in values]

print(columns[:4],'\n')

for i in reformated_values[:4]:
	print(i[:4],'\n')

print(len(columns),'\n')
# # print(columns,'\n')

print(len(reformated_values),'\n')
# # print(reformated_values[19:22])


# # Construct a dataframe from genes of interest, 
# # with genes as rows and samples as columns.

Df = pd.DataFrame(reformated_values, columns=columns)


print(Df.iloc[:100,:10])


# # Save the formatted gene set and control set dataframes as CSV files.

# Df.to_csv('../../../Data/RNA_Data/TCGA_mRNA_TPM.csv', index = False)


