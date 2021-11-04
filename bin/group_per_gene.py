#!/usr/bin/env python
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, required=True, help='A per target coverage bedfile that is a outputfile from gvcf2bed2.py')
parser.add_argument('--output', type=str, required=True, help='output filename with calculations per gene. (.bed)')
args = parser.parse_args()

import_file_per_target=args.input
output_file_per_gene=args.output

'''
read in targetfile: *.merged.variant.calls.per_gene.bed
'''
df = pd.read_csv(import_file_per_target,sep='\t')
print(df.head())

'''
remove '%' , and cast to float
'''
df['GQ_low']    = df['GQ_low'].replace({'%': ''}, regex=True)
df['GQ_low']    = df['GQ_low'].astype(float)
df['GQ_medium'] = df['GQ_medium'].replace({'%': ''}, regex=True)
df['GQ_medium'] = df['GQ_medium'].astype(float)
df['GQ_high']   = df['GQ_high'].replace({'%': ''}, regex=True)
df['GQ_high']   = df['GQ_high'].astype(float)
df.head()

'''
get gene start en stop
'''
gene_start = df.groupby('gene')['start'].min()
gene_stop = df.groupby('gene')['stop'].max()
df['gene_start'] = df['gene'].map(gene_start)
df['gene_stop'] = df['gene'].map(gene_stop)
print(df.head())

'''
Group per gene, and ignore gene_start gene_stop
'''

final_df=df.groupby(['gene','chr'],as_index=False).agg({
'gene_start': np.min,
'gene_stop': np.max,
'GQ_low': np.mean,        # is this correct?
'GQ_medium': np.mean,     # is this correct?
'GQ_high': np.mean,       # is this correct?
'avgDp': np.mean,
'percentage DP<20': np.mean,
'percentage DP>20': np.mean,
'summedDP': np.sum,
'targetSize': np.sum,
'numberBasesDPInCategoryLow': np.sum,
'numberBasesZeroCoverage': np.sum
})

'''
roundup floats
'''
final_df['GQ_low']=final_df['GQ_low'].round(2)
final_df['GQ_medium']=final_df['GQ_medium'].round(2)
final_df['GQ_high']=final_df['GQ_high'].round(2)
final_df['avgDp']=final_df['summedDP']/final_df['targetSize'] #recalulate avgDp
final_df['avgDp']=final_df['avgDp'].round(2)
final_df['percentage DP<20']=final_df['percentage DP<20'].round(2)
final_df['percentage DP>20']=final_df['percentage DP>20'].round(2)
print(final_df.head())

'''
output csv
'''
final_df.to_csv(output_file_per_gene,index=False,sep='\t')
print('output file: ' + output_file_per_gene )
