import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import statsmodels.stats.multitest as smm
from scipy.stats import ttest_ind
from scipy.stats import mannwhitneyu
import numpy as np

sns.set_theme(style="whitegrid")

pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)


"""INITIAL DATASET CONFIGURATION
"""

# Loading the datasets with Pandas. 
cntrl_RNA_df = pd.read_csv('../../Data/RNA_Data/Control_Genes_RNA_df.csv')
GoI_RNA_df = pd.read_csv('../../Data/RNA_Data/Genes_of_Interest_RNA_df.csv')
cntrl_CNV_df = pd.read_csv('../../Data/CNV_Data/Control_Genes_CNV_df.csv')
GoI_CNV_df = pd.read_csv('../../Data/CNV_Data/Genes_of_Interest_CNV_df.csv')
codes = pd.read_csv('../../Data/Codes.csv')
hist_df = pd.read_csv('../../Data/Cancer_Groups.csv')



codes['TSS Code'] = codes['TSS Code'].astype(str)
codes['TSS Code'] = codes['TSS Code'].apply(lambda x: '0' + x if len(x) == 1 else x)
codes['TSS Code'] = codes['TSS Code'].replace('nan', 'NA')


TSS_mapping_dict = pd.Series(codes['Study Name'].values, index=codes['TSS Code']).to_dict()
Hist_mapping_dict = pd.Series(hist_df['Type'].values, index=hist_df['Name']).to_dict()
TSS_abbr_mapping_dict = pd.Series(hist_df['Abbrvs'].values, index=hist_df['Name']).to_dict()



"""PART 1: INITIAL DISTRIBUTION ANALYSIS

PURPOSE: Add broad metric variables (Mean, Median, Standard deviation) 
to the data and analyse the general distribution of the data based on 
group type.

"""


# Create the columns for analysis metrics in the
# set and control dataframes.

Stats_SoI_RNA_df = GoI_RNA_df.iloc[:,:]
Stats_Cntrl_RNA_df = cntrl_RNA_df.iloc[:,:]
print(Stats_SoI_RNA_df.iloc[:5,-5:])

Stats_SoI_CNV_df = GoI_CNV_df.iloc[:,:]
Stats_Cntrl_CNV_df = cntrl_CNV_df.iloc[:,:]
print(Stats_SoI_CNV_df.iloc[:5,-5:])


Stats_SoI_CNV_df['Mean'] = Stats_SoI_CNV_df.iloc[:, 1:].mean(axis=1)
Stats_SoI_CNV_df['Median'] = Stats_SoI_CNV_df.iloc[:, 1:-1].median(axis=1)
Stats_SoI_CNV_df['Std'] = Stats_SoI_CNV_df.iloc[:, 1:-2].std(axis=1)

Stats_Cntrl_CNV_df['Mean'] = Stats_Cntrl_CNV_df.iloc[:, 1:].mean(axis=1)
Stats_Cntrl_CNV_df['Median'] = Stats_Cntrl_CNV_df.iloc[:, 1:-1].median(axis=1)
Stats_Cntrl_CNV_df['Std'] = Stats_Cntrl_CNV_df.iloc[:, 1:-2].std(axis=1)


Stats_SoI_RNA_df['Mean'] = Stats_SoI_RNA_df.iloc[:, 1:].mean(axis=1)
Stats_SoI_RNA_df['Median'] = Stats_SoI_RNA_df.iloc[:, 1:-1].median(axis=1)
Stats_SoI_RNA_df['Std'] = Stats_SoI_RNA_df.iloc[:, 1:-2].std(axis=1)

Stats_Cntrl_RNA_df['Mean'] = Stats_Cntrl_RNA_df.iloc[:, 1:].mean(axis=1)
Stats_Cntrl_RNA_df['Median'] = Stats_Cntrl_RNA_df.iloc[:, 1:-1].median(axis=1)
Stats_Cntrl_RNA_df['Std'] = Stats_Cntrl_RNA_df.iloc[:, 1:-2].std(axis=1)
 

# Create column defining the groups.

Stats_SoI_RNA_df['Group'] = 'Set'
Stats_Cntrl_RNA_df['Group'] = 'Control'

Stats_SoI_CNV_df['Group'] = 'Set'
Stats_Cntrl_CNV_df['Group'] = 'Control'


# Reduce the dataframes to the necessary columns.

reduced_SoI_RNA = Stats_SoI_RNA_df[['Sample', 'Group', 'Mean','Median', 'Std']]
reduced_cntrl_RNA = Stats_Cntrl_RNA_df[['Sample', 'Group', 'Mean','Median', 'Std']]

reduced_SoI_CNV = Stats_SoI_CNV_df[['Sample', 'Group', 'Mean','Median', 'Std']]
reduced_cntrl_CNV = Stats_Cntrl_CNV_df[['Sample', 'Group', 'Mean','Median', 'Std']]


## FULL SAMPLE CREATION


# Merge the reduced dataframes for further analysis in R 

reduced_RNA_full = pd.merge(reduced_SoI_RNA, reduced_cntrl_RNA, how = 'inner', on= 'Sample', suffixes=('_set', '_cntrl'))
reduced_CNV_full = pd.merge(reduced_SoI_CNV, reduced_cntrl_CNV, how = 'inner', on= 'Sample', suffixes=('_set', '_cntrl'))


# Add the tissue source type to the dataframe.

reduced_RNA_full['TSS'] = reduced_RNA_full.Sample.apply(lambda x: x.split('-')[1])
reduced_RNA_full['TSS'] = reduced_RNA_full['TSS'].map(TSS_mapping_dict)

reduced_CNV_full['TSS'] = reduced_CNV_full.Sample.apply(lambda x: x.split('-')[1])
reduced_CNV_full['TSS'] = reduced_CNV_full['TSS'].map(TSS_mapping_dict)


# Add the Histology Group to the dataframe.

reduced_RNA_full['Hist_group'] = reduced_RNA_full['TSS'].map(Hist_mapping_dict)
reduced_RNA_full['TSS_Abbrvs'] = reduced_RNA_full['TSS'].map(TSS_abbr_mapping_dict)

reduced_CNV_full['Hist_group'] = reduced_CNV_full['TSS'].map(Hist_mapping_dict)
reduced_CNV_full['TSS_Abbrvs'] = reduced_CNV_full['TSS'].map(TSS_abbr_mapping_dict)

print(reduced_CNV_full.iloc[:5,:5])

print(reduced_RNA_full.iloc[:5,:5])

# Save the datasets.

# reduced_RNA_full.to_csv('../../Data/RNA/Full_RNA_Sample_metrics.csv', index=False)
# reduced_CNV_full.to_csv('../../Data/CNV/Full_CNV_Sample_metrics.csv', index=False)


