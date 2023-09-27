#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 12:52:57 2023

@author: francois
"""
# =============================================================================
# #### Pairplots for selection coefficient
# =============================================================================
import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.stats as stats
import scipy
import seaborn as sns
from collections import Counter
import re
from Bio import SeqIO
from subprocess import call
import math
import matplotlib.patches as mpatches
import sys
import statistics
from scipy.stats import spearmanr
from scipy.stats import pearsonr
#%%
runcell(0)
aa_list = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", 
                "A", "V", "I", "L", "M", "F", "Y", "W", "*"]
seq ="DWQKSLTLIVALTLSRGIGLKNDLPWKLKSDMMFFSRVTSGLLVTRSTGQMNVVLMGRKTWESLPAHSRPLKNRINVVISRQEVLDLGGGAYHARSLDDALALLSQIYDSTSKIQLNRVFVIGGGELCKAAMEHSRLNRIIATVIHNEVDCDVFFPIDFRSSQSCLPWRKQDHSVLEAWVGSKVPQGKINENGFIYEFEMWIRD"
positions = list(range(2,206))

# Run cell to load condition data
condition = "MTX20" #"DMSO" for DMSO, "MTX4" for IC75 and "MTX20" for IC90
longform = pd.read_csv("../Longform_dataframes/total_"+condition+"_longform_CoS_t15.csv")


# Create silent dataframe as reference
longform_silent = longform.copy()
longform_silent = longform_silent[longform_silent['nature']=='silent']
longform_silent = longform_silent.reset_index()
cols_to_keep = [condition+'-1_CoS', condition+'-2_CoS', condition+'-3_CoS', 'index']
cols_to_keep2 = [condition+'-1_CoS', condition+'-2_CoS', condition+'-3_CoS']

longform_silent_melt = longform_silent.copy()
longform_silent_melt = longform_silent_melt[cols_to_keep]
longform_silent_melt = longform_silent_melt.melt(id_vars=['index'], value_vars = cols_to_keep, ignore_index = False, var_name="Sample", value_name=condition+'_CoS')

longform['p-values'] = 0.00
longform['Welch_t-test'] = 0.00

silent_array = longform_silent_melt[condition+'_CoS']

longform_aa = pd.DataFrame(columns=['aa_WT','position','aa_mut','p-values','Welch_t-test','codons_num','median_CoS','fragment',condition+'-1_CoS', condition+'-2_CoS', condition+'-3_CoS'])
index = 0

# Run code to aggregated positions by AA and run statistical tests

for posit in positions:
    mutant_df_master = longform[longform['position'] == posit]
    fragment = longform.at[index,"fragment"]
    for aa in aa_list:
        aa_df = mutant_df_master[mutant_df_master['aa_mut'] == aa]
        if len(aa_df) > 0:
            aa_wt = seq[posit-2]
            mutant_df_array = aa_df[cols_to_keep2]
            mutant_df_median = aa_df['median_CoS'].median()
            mutant_df_median1 = aa_df[condition+'-1_CoS'].median()
            mutant_df_median2 = aa_df[condition+'-2_CoS'].median()
            mutant_df_median3 = aa_df[condition+'-3_CoS'].median()
            mutant_array = mutant_df_array.melt(value_vars = cols_to_keep2, ignore_index = True)
            mutant_array = mutant_array['value']

            t, p = stats.ttest_ind(silent_array, mutant_array, equal_var = False)
            longform_aa.at[index,'aa_WT'] = aa_wt
            longform_aa.at[index,'position'] = posit
            longform_aa.at[index,'aa_mut'] = aa
            longform_aa.at[index,'p-values'] = p
            longform_aa.at[index,'Welch_t-test'] = t
            longform_aa.at[index,'codons_num'] = len(mutant_array)/3
            longform_aa.at[index,'median_CoS'] = mutant_df_median
            longform_aa.at[index,'fragment'] = fragment
            longform_aa.at[index,condition+'-1_CoS'] = mutant_df_median1
            longform_aa.at[index,condition+'-2_CoS'] = mutant_df_median2
            longform_aa.at[index,condition+'-3_CoS'] = mutant_df_median3
            index = index + 1
        else:
            pass
longform_aa = longform_aa.astype({'aa_WT': str, 
                                  'position': 'int64',
                                  'aa_mut': 'object',
                                  'p-values': 'float64',
                                  'Welch_t-test': 'float64',
                                  'codons_num': 'float64',
                                  'median_CoS': 'float64'}, copy = True)

# Generate plots to compare different replicatesi in given conditions (Figure S1ABC)
if condition == "DMSO":
    IC = "DMSO"
    supp = 'Supp3A'
    filename = "./FigureS3A"
    rou = 2
elif condition == "MTX4":
    IC = 'IC75'
    supp = 'Supp3B'
    filename = "./FigureS3B"
    rou = 2
elif condition == "MTX20":
    IC = 'IC90'
    supp = 'Supp3C'
    filename = "./FigureS3C"
    rou = 2

def corrfunc(x, y, **kws):
    (r, p) = spearmanr(x, y)
    p = round(p,rou)
    ax = plt.gca()
    if r >= 0.999:
        pass
    else:
        ax.annotate("$\itρ$ = {:.2f} ".format(r),
                    xy=(0.05, .9), xycoords=ax.transAxes)
        ax.annotate("p < 1e-10",
                    xy=(0.05, .82), xycoords=ax.transAxes)


cols = [condition+'-1_CoS', condition+'-2_CoS',condition+'-3_CoS']
df = longform_aa.copy()

df = df.loc[:, cols]
renamer = {condition+'-1_CoS' : 'Replicate 1',condition+'-2_CoS':'Replicate 2', condition+'-3_CoS' :'Replicate 3'}
df.rename(columns = renamer, inplace = True)
df = df.astype('float64')
plot = sns.pairplot(data = df,
                    corner = True, kind = 'reg', height=3.333)
plot.map(corrfunc)
plot.axes[0,0].set_xlim(-1.5,2)
plot.axes[0,0].set_ylim(-1.5,2)
plot.axes[1,0].set_xlim(-1.5,2)
plot.axes[1,0].set_ylim(-1.5,2)
plot.axes[1,1].set_xlim(-1.5,2)
plot.axes[1,1].set_ylim(-1.5,2)
plot.axes[2,0].set_xlim(-1.5,2)
plot.axes[2,0].set_ylim(-1.5,2)
plot.axes[2,1].set_xlim(-1.5,2)
plot.axes[2,1].set_ylim(-1.5,2)
plot.axes[2,2].set_xlim(-1.5,2)
plot.axes[2,2].set_ylim(-1.5,2)

plt.savefig(filename, dpi=200, bbox_inches="tight")
#%%
# Figure S8
aa_list = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", 
                "A", "V", "I", "L", "M", "F", "Y", "W", "*"]
seq ="DWQKSLTLIVALTLSRGIGLKNDLPWKLKSDMMFFSRVTSGLLVTRSTGQMNVVLMGRKTWESLPAHSRPLKNRINVVISRQEVLDLGGGAYHARSLDDALALLSQIYDSTSKIQLNRVFVIGGGELCKAAMEHSRLNRIIATVIHNEVDCDVFFPIDFRSSQSCLPWRKQDHSVLEAWVGSKVPQGKINENGFIYEFEMWIRD"
positions = list(range(2,206))

# Run cell to load condition data
longform = pd.read_csv("../Longform_dataframes/total_DMSO_longform_CoS_t15.csv")

# Create silent dataframe as reference
longform_silent = longform.copy()
longform_silent = longform_silent[longform_silent['nature']=='silent']
longform_silent = longform_silent.reset_index()
condition = "DMSO"
cols_to_keep = [condition+'-1_CoS', condition+'-2_CoS', condition+'-3_CoS', 'index']
cols_to_keep2 = [condition+'-1_CoS', condition+'-2_CoS', condition+'-3_CoS']

longform_silent_melt = longform_silent.copy()
longform_silent_melt = longform_silent_melt[cols_to_keep]
longform_silent_melt = longform_silent_melt.melt(id_vars=['index'], value_vars = cols_to_keep, ignore_index = False, var_name="Sample", value_name=condition+'_CoS')

longform['p-values'] = 0.00
longform['Welch_t-test'] = 0.00

silent_array = longform_silent_melt[condition+'_CoS']

longform_aa = pd.DataFrame(columns=['aa_WT','position','aa_mut','p-values','Welch_t-test','codons_num','median_CoS','fragment',condition+'-1_CoS', condition+'-2_CoS', condition+'-3_CoS'])
index = 0

# Run code to aggregated positions by AA and run statistical tests

for posit in positions:
    mutant_df_master = longform[longform['position'] == posit]
    fragment = longform.at[index,"fragment"]
    for aa in aa_list:
        aa_df = mutant_df_master[mutant_df_master['aa_mut'] == aa]
        if len(aa_df) > 0:
            aa_wt = seq[posit-2]
            mutant_df_array = aa_df[cols_to_keep2]
            mutant_df_median = aa_df['median_CoS'].median()
            mutant_df_median1 = aa_df[condition+'-1_CoS'].median()
            mutant_df_median2 = aa_df[condition+'-2_CoS'].median()
            mutant_df_median3 = aa_df[condition+'-3_CoS'].median()
            mutant_array = mutant_df_array.melt(value_vars = cols_to_keep2, ignore_index = True)
            mutant_array = mutant_array['value']

            t, p = stats.ttest_ind(silent_array, mutant_array, equal_var = False)
            longform_aa.at[index,'aa_WT'] = aa_wt
            longform_aa.at[index,'position'] = posit
            longform_aa.at[index,'aa_mut'] = aa
            longform_aa.at[index,'p-values'] = p
            longform_aa.at[index,'Welch_t-test'] = t
            longform_aa.at[index,'codons_num'] = len(mutant_array)/3
            longform_aa.at[index,'median_CoS'] = mutant_df_median
            longform_aa.at[index,'fragment'] = fragment
            longform_aa.at[index,condition+'-1_CoS'] = mutant_df_median1
            longform_aa.at[index,condition+'-2_CoS'] = mutant_df_median2
            longform_aa.at[index,condition+'-3_CoS'] = mutant_df_median3
            index = index + 1
        else:
            pass
longform_aa = longform_aa.astype({'aa_WT': str, 
                                  'position': 'int64',
                                  'aa_mut': 'object',
                                  'p-values': 'float64',
                                  'Welch_t-test': 'float64',
                                  'codons_num': 'float64',
                                  'median_CoS': 'float64'}, copy = True)

stop = longform[longform['nature']=='stop']
plt.figure(figsize = (10,10))
sns.scatterplot(stop, x = "position", y = "median_CoS", color = 'black')
plt.ylabel('Selection coefficient - DMSO', fontsize = 20)
plt.xlabel('Position', fontsize = 20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.savefig('./FigureS8.png', dpi=200, bbox_inches="tight")
