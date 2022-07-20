#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import datasets, linear_model
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats import weightstats as stests
from scipy.stats import ttest_ind

def analyze_for_drug(df, drug_name):
    df.hist(column=['Dose', 'GFP_nuc_to_cyto_ratio'])
    df.plot.scatter(x='Dose', y='GFP_nuc_to_cyto_ratio')
    print(" --------- {:} --------".format(drug_name))
    print("Ratio by Column: ", df.groupby('column').GFP_nuc_to_cyto_ratio.mean())
    print("Ratio by Dose: ", df.groupby('Dose').GFP_nuc_to_cyto_ratio.mean())

def hist_per_dose(df, title, subfolder):
    doses = df.Dose.unique()
    df_per_dose = []
    for dose in doses:
        df_per_dose.append(df[ df.Dose == dose ])
    
    for dft in df_per_dose:
        dose = dft.Dose.unique()[0]
        well_column = dft.column.unique()[0]
        print("Dose = ", dose)
        print("Column = ", well_column)
        ax = dft.hist(column=['GFP_nuc_to_cyto_ratio'], figsize=(20, 10))
        fig = ax[0,0].get_figure()
        fig_filename = "./output_pics/{:}/hist_per_dose_{:}_{:02}_{:}.png".format(subfolder, title, well_column, dose)
        fig.savefig(fig_filename)
        #print("GFP_nuc_to_cyto_ratio", dft.GFP_nuc_to_cyto_ratio)
        #print(np.histogram(dft.GFP_nuc_to_cyto_ratio))

def threshold_GFP_ratio(df):
    dft = df[df.GFP_nuc_to_cyto_ratio>0]
    return dft

def plot_scatter_and_curves(df, titlestr):
   
    plt.rcParams["figure.figsize"] = (20,10)
    
    x = df.Dose.values
    y = df.GFP_nuc_to_cyto_ratio.values    
    plt.plot(x, y, '.')
    
    dft = df.groupby('Dose')['GFP_nuc_to_cyto_ratio'].mean()    
    x = dft.index.values
    avgs = dft.values
    ax=plt.plot(x, avgs)
    plt.xlabel('Dose', fontsize=18)
    plt.ylabel('Response Ratio', fontsize=18)
    plt.title(titlestr, fontsize=18)
    plt.show()
    
def get_correlation(df, methodname = 'pearson'):
    # {‘pearson’, ‘kendall’, ‘spearman’}
    corr = df.corr(method = methodname)    
    print (corr.iloc[5])
    #sns.heatmap(corr)

def compare_druggroups(drug1, drug2):
    ttest,pval = ttest_ind(drug1.GFP_nuc_to_cyto_ratio.values, drug2.GFP_nuc_to_cyto_ratio.values, equal_var=False)
    print ("t-test-ind: pval = ", pval)

def analyze_per_drug(df_wortmannin, df_ly):
    plot_scatter_and_curves(df_wortmannin, 'Wortmannin')
    plot_scatter_and_curves(df_ly, 'LY294.002')
    
    get_correlation(df_wortmannin, 'spearman')
    get_correlation(df_ly, 'spearman')
    
    compare_druggroups(df_wortmannin, df_ly)    
    
# ------------------------------------------------------
    
dfbase = pd.DataFrame.from_csv('gfp_cells_complete.csv')

df_wortmannin = dfbase[(dfbase.row >= 'A') & (dfbase.row <= 'D') & (dfbase.column >= 3) & (dfbase.column <= 11)]
df_ly = dfbase[(dfbase.row >= 'E') & (dfbase.row <= 'H') & (dfbase.column >= 3) & (dfbase.column <= 11)]
#hist_per_dose(df_wortmannin, 'wortmannin', 'no_threshold')
#hist_per_dose(df_ly, 'ly', 'no_threshold')

# remove cells with no transport in nucleus *
df_wortmannin = threshold_GFP_ratio(df_wortmannin)
df_ly = threshold_GFP_ratio(df_ly)

# analyze -------------------------------------------
analyze_per_drug(df_wortmannin, df_ly)

# include the wells with no drug added -------------------------------
df_no_drug = dfbase[ dfbase.Dose == 0 ]
df_no_drug = threshold_GFP_ratio(df_no_drug)
df_wortmannin = df_wortmannin.append(df_no_drug)
df_ly = df_ly.append(df_no_drug)

analyze_per_drug(df_wortmannin, df_ly)

#analyze_for_drug(df_wortmannin, 'Wortmannin')
#analyze_for_drug(df_ly, 'LY294.002')

