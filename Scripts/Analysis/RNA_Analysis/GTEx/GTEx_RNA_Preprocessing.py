

'''

This Script is a python script which serves to clean and preprocess the
GTEx RNA data that contains the healthy samples that are to be used in 
the Part one of the Analysis. 

The first part of the cleaning will be to remove the samples (subjects)
which do not meet the criteria of being healthy. 
This is based on the DTHHRDY (Death Hardy Scale) measure from GTEx. 
Which is contained in the meta file that was also acquired. This
measure is a value ranging from 0-4 and pertains to the state of the 
subject at the time of or before death. Values 0, 1, 2, and 3 are suitable for 
the analysis as they refer to cases involving sudden deaths (<24hrs) or ventilator cases that 
were not due to an underlying disease. A value of 4 is usually used
for deaths in which cancer was involved. 


'''



import pandas as pd
import dask.dataframe as dd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import statsmodels.stats.multitest as smm


pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)


# Loading of the full RNA GTEx dataset.
GTEx_Samp_df = dd.read_csv(
	'../../../Data/RNA_Data/GTEx_RNA_Healthy/GTEx_RNA_Full.csv',
    storage_options={'anon': True, 'use_ssl': False},
    sample = 100000000
    )



# Loading of the full RNA GTEx Annotation dataset.
GTEx_meta_df = pd.read_csv('../../../Data/RNA_Data/GTEx_RNA_Healthy/GTEx_RNA_Subject_meta.csv')

# print(GTEx_Samp_df.shape[1])
# The Curent number of samples in the GTEx Dataset stands at 17385

# Get the SUBJIDS that correspond to the samples to be used
Healthy_Subjects = GTEx_meta_df[GTEx_meta_df.DTHHRDY.isin([ 1, 2, 3, 0])].SUBJID.unique()

# Get the column names from the full RNA GTEX dataset.
cols = list(GTEx_Samp_df.columns[3:])

# Get the columsn that should be dropped from the dataset.
cols_to_drop = [subj for subj in cols if subj[:10] not in Healthy_Subjects]

# Drop the columns.
GTEx_Samp_df = GTEx_Samp_df.drop(columns = cols_to_drop)
GTEx_Samp_df = GTEx_Samp_df.drop(columns = ['Unnamed: 0'])

# Save the new dataset
# GTEx_Samp_df.to_csv('../../../Data/RNA_Data/GTEx_RNA_Healthy/GTEx_RNA_Healthy*.csv',
#     index = False)



