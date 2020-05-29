import numpy as np 
import pandas as pd
import sys
import gzip
import urllib.request
from sklearn.metrics import roc_curve    # for first part
from sklearn.metrics import roc_auc_score # for first part 
from matplotlib import pyplot as plt
sys.path.append('..')
from utils import tools as tl
import bioinfokit  # for volcano plots
from bioinfokit import visuz

# load .vcf file in GHr38 format
df38 = tl.read_vcf('..\\data\\VCF_phenotypes_SNPs\\diabetes.vcf')

# drop duplicates 
df38.drop_duplicates(inplace=True)

# process
df38 = df38[~df38['CHROM'].str.contains('CHR')]
df38 = df38[~df38['ALT'].str.contains(',')]

# process dataframe for Liftover conversion (GHr38->GHr37)
df37 = df38.drop(columns=['INFO', 'FILTER', 'QUAL', 'ALT', 'REF', 'ID'])
x = df37['POS'].to_numpy()
y = x + 1
df37['POS+1'] = y
df37['CHROM'] = 'chr' + df37['CHROM'].astype(str)

# *** FILE GOES TO LIFTOVER AT THIS POINT *** #
tl.to_DeepSEAbed(df37, 'temp\\Gr38_.bed')
df37_new = tl.load_bed('temp\\GR38_converted.bed')

df37_new['ref']=df38['REF'].to_numpy()
df37_new['alt']=df38['ALT'].to_numpy()
df37_new['chromEnd'] = df38['ID'].to_numpy()
df37_new.rename(columns = {"chromEnd" : "ID"}, inplace=True)



tl.to_DeepSEAbed(df37_new, '..\\data\\VCF_test_data\\diabetes_variants.vcf')
results = tl.load_out_bed('..\\data\\VCF_test_output\\diabetes\\infile.vcf.out.logfoldchange')
e_values = x = tl.load_out_bed('..\\data\\VCF_test_output\\diabetes\\infile.vcf.out.evalue')
  

# Volcano Plot for a particular SNP

snp_df = pd.DataFrame(columns = ['Names', 'log2FC', 'E-value'])
snp_e_values = e_values.drop(columns = ['chrom', 'pos', 'name', 'ref', 'alt']).iloc[204].to_numpy()
snp_logfold = results.drop(columns = ['chrom', 'pos', 'name', 'ref', 'alt']).iloc[204].to_numpy()
snp_ids = pd.DataFrame(columns=e_values.columns).drop(columns=['chrom', 'pos', 'name', 'ref', 'alt'])
snp_ids = snp_ids.columns

snp_df['log2FC'] = snp_logfold
snp_df['E-value'] = snp_e_values
snp_df['Names'] = snp_ids

snp_df
tf_labels = tl.extract_features(tf=True)
els = cols = [x for x in tf_labels if 'USF1' in x]
#els_2 = cols = [x for x in tf_labels if 'CEBPB'  in x]
visuz.gene_exp.volcano(d=snp_df, lfc="log2FC", pv="E-value", geneid="Names", genenames=({els[0]:'USF1'}), gstyle=2)


results = tl.load_out_bed('..\\data\\VCF_test_output\\diabetes\\infile.vcf.out.logfoldchange')
results = results.drop(columns = ['chrom', 'pos', 'name', 'ref', 'alt'])

cell_cols = [col for col in results.columns if 'K562' in col]
result_cell = results[cell_cols]
feature_cols = [col for col in result_cell.columns if 'H3' in col or 'Nrf' in col]

x = result_cell[feature_cols]

df = x
visuz.stat.corr_mat(df=df)
