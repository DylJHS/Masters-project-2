import pandas as pd
import dask.dataframe as dd
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.dataset as ds
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import statsmodels.stats.multitest as smm

sns.set_theme(style="whitegrid")

pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

cntrl_CNV_df = pd.read_csv('../../Data/CNV/Control_Genes_CNV_df.csv')
GoI_CNV_df = pd.read_csv('../../Data/CNV/Genes_of_Interest_CNV_df.csv')

reduced_RNA_metrics = pd.read_csv('../../Data/RNA/Full_RNA_Sample_metrics.csv')
reduced_CNV_metrics = pd.read_csv('../../Data/CNV/Full_CNV_Sample_metrics.csv')

print(cntrl_CNV_df.shape)
print(GoI_CNV_df.shape)

print(reduced_RNA_metrics.shape)
print(reduced_CNV_metrics.shape)

