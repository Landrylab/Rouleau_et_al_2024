#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 18:16:58 2023

@author: francois
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math
import scipy.stats as stats
from scipy.stats import spearmanr
from matplotlib.ticker import FormatStrFormatter

#%%
# Set condition you want to analyse

# =============================================================================
condition = 'DMSO' # Or MTX4 (for IC75) Or MTX20 (for IC90)
# =============================================================================
data = pd.read_csv('../Longform_dataframes/total_'+condition+'_longform_CoS_t15.csv')
back = pd.read_excel("./CoS_dgr_replicates_error.xlsx")
to_fill = back.copy()

#%%
# Create smaller dataframe with only mutants of interest and fill
for row in to_fill.index:
    position = to_fill.at[row,'position']
    codon_mut = to_fill.at[row,'codon']
    if codon_mut == 'TAA':
        codon_mut = 'TAG'
        temp_df = data[data['position']==position]
        temp_df = temp_df[temp_df['codon_mut'] == codon_mut].reset_index()
        median = temp_df.at[0,'median_CoS']
        rep1 = temp_df.at[0,condition+'-1_CoS']
        rep2 = temp_df.at[0,condition+'-2_CoS']
        rep3 = temp_df.at[0,condition+'-3_CoS']
        
        to_fill.at[row,'median_'+condition] = median
        to_fill.at[row,'rep1_'+condition] = rep1
        to_fill.at[row,'rep2_'+condition] = rep2
        to_fill.at[row,'rep3_'+condition] = rep3

# Save dataframe with errorbar information
to_fill.to_csv("./CoS_dgr_replicates_error_temp.csv", sep = ',')
    
#%%
# This will take a while to reget DGR for all replicates
# Might actually have to fill excel by hand
back = pd.read_csv("./CoS_dgr_replicates_error_full.csv", index_col = 0)
data = back.copy()
data[['MTX20_max','MTX20_min','MTX4_max','MTX4_min','DMSO_max','DMSO_min','MTX_dgr_max','MTX_dgr_min','DMSO_dgr_max','DMSO_dgr_min']] = np.nan
data = data.reset_index(drop = True)
#%%
# Calculate errorbars
medians = ['median_DMSO', 'median_MTX4', 'median_MTX20', 'median_dgr-MTX', 'median_dgr-DMSO']


for row in data.index:
    data_mutant = data.iloc[row,]
    DMSO_median = data_mutant.at['median_DMSO']
    array_CoS_DMSO = [data_mutant.at['rep1_DMSO'],data_mutant.at['rep2_DMSO'],data_mutant.at['rep3_DMSO']]
    
    MTX4_median = data_mutant.at['median_MTX4']
    array_CoS_MTX4 = [data_mutant.at['rep1_MTX4'],data_mutant.at['rep2_MTX4'],data_mutant.at['rep3_MTX4']]
    
    MTX20_median = data_mutant.at['median_MTX20']
    array_CoS_MTX20 = [data_mutant.at['rep1_MTX20'],data_mutant.at['rep2_MTX20'],data_mutant.at['rep3_MTX20']]
    
    MTXdgr_median = data_mutant.at['median_dgr_MTX']
    array_dgr_MTX = [data_mutant.at['dgr1_MTX'],data_mutant.at['dgr2_MTX'],data_mutant.at['dgr3_MTX']]
    
    DMSOdgr_median = data_mutant.at['median_dgr_DMSO']
    array_dgr_DMSO = [data_mutant.at['dgr1_DMSO'],data_mutant.at['dgr2_DMSO'],data_mutant.at['dgr3_DMSO']]

    array_array = [array_CoS_DMSO,array_CoS_MTX4,array_CoS_MTX20,array_dgr_MTX,array_dgr_DMSO]
    
    for array in array_array:
        if array == array_CoS_DMSO:
            median = DMSO_median
            for value in array:
                if value > median:
                    up = value - median
                    data.at[row,'DMSO_max'] = up
                if value < median:
                    down = median - value
                    data.at[row,'DMSO_min'] = down
        if array == array_CoS_MTX4:
            median = MTX4_median
            for value in array:
                if value > median:
                    up = value - median
                    data.at[row,'MTX4_max'] = up
                if value < median:
                    down = median - value
                    data.at[row,'MTX4_min'] = down
        if array == array_CoS_MTX20:
            median = MTX20_median
            for value in array:
                if value > median:
                    up = value - median
                    data.at[row,'MTX20_max'] = up
                if value < median:
                    down = median - value
                    data.at[row,'MTX20_min'] = down
        if array == array_dgr_MTX:
            median = MTXdgr_median
            for value in array:
                if value > median:
                    up = value - median
                    data.at[row,'MTX_dgr_max'] = up
                if value < median:
                    down = median - value
                    data.at[row,'MTX_dgr_min'] = down
        if array == array_dgr_DMSO:
            median = DMSOdgr_median
            for value in array:
                if value > median:
                    up = value - median
                    data.at[row,'DMSO_dgr_max'] = up
                if value < median:
                    down = median - value
                    data.at[row,'DMSO_dgr_min'] = down
  
  
data = data.fillna(0)
  

#%%
# Panel A
# IC75 VS IC90
# Set condition
# =============================================================================
condition = 'MTX' 
# =============================================================================
# Make backup
df = data.copy()
# Drop NA
df_ok = df.dropna()

# Set data/parameters of interest depending on which condition you use
if condition == "MTX":
    condition_num = 'MTX20'
    to_keep = ['mutant','median_MTX20', 'MTX20_max', 'MTX20_min', 'median_dgr_MTX', 'MTX_dgr_max', 'MTX_dgr_min']
    #df = data[to_keep]
    height = -0.7
    pal = "coolwarm"
    ex = 0.7
    IC = 'IC90'
    empty = 0.0459
    sc = 0.196
    wt = "#82A6FB"
    pp = 14
elif condition == "DMSO":
    to_keep = ['mutant','median_DMSO', 'DMSO_max', 'DMSO_min', 'median_dgr_DMSO', 'DMSO_dgr_max', 'DMSO_dgr_min']
    #df = data[to_keep]
    height = 0.2
    pal = sns.diverging_palette(135, 30, l=65, center="dark", as_cmap=True)
    ex = -0.9
    empty = 0.3075
    sc = 0.289
    wt = "#151515"
    pp = 5

#Drop control values from dataframe
df_ex = df_ok[df_ok["mutant"] != "PjDHFR"]
df_ex = df_ex[df_ex["mutant"] != "Empty"]


# Generate figure
fig, ax = plt.subplots(figsize = (6,6))

# Create scatterplot with color palette
out = ['G58Q', 'D110W', 'S111E', 'I196W', 'I196F']
g = sns.scatterplot(data=df_ex, x='median_MTX20', y='median_MTX4', hue='median_MTX20', palette=pal, vmin=-0.5, vmax=1.5, legend=False, zorder=10, ax=ax,linewidth = 0)

# Define and plot error bars
xerror = [df_ex['MTX4_min'].to_list(),df_ex['MTX4_max'].to_list()]
yerror = [df_ex['MTX20_min'].to_list(),df_ex['MTX20_max'].to_list()]
g.errorbar(x=df_ex['median_MTX20'], y=df_ex['median_MTX4'], xerr = xerror, yerr = yerror, color = 'black', fmt = 'o')
g.axvline

# Highlight the 'wt' data point
df_wt = df_ok[df_ok["mutant"] == "PjDHFR"]
g1 = sns.scatterplot(df_wt, x='median_MTX20', y='median_MTX4', color=wt, marker='X', s=100, legend=False, zorder=10,ax=ax,linewidth = 0)

# Create regplot
g2 = sns.regplot(data=df_ex, x='median_MTX20', y='median_MTX4', scatter=False, color='lightgrey',ax=ax)

# Plot colorbar
norm = plt.Normalize(-0.75, 1.5)
sm = plt.cm.ScalarMappable(cmap=pal,norm = norm)
cax = fig.add_axes([ax.get_position().x1+0.1, ax.get_position().y0+0.04, 0.02, ax.get_position().height])
ax.figure.colorbar(sm, cax=cax)


df_out = df.loc[df['mutant'].isin(out)]
g3 = sns.scatterplot(data=df_out, x='median_MTX20', y='median_MTX4', marker='o',ax=ax, zorder=9, s= 100, linewidth = 0, color = "green")



# Add other plot elements
coef, p = spearmanr(df_ok['median_MTX20'], df_ok['median_MTX4'])
coef = round(coef, 2)
p = round(p, 22)
g1.text(s="$\itρ$ = {}\np = {}".format(coef, p), x=ex, y=height, fontsize=20)
g1.set_xlabel("Selection coefficient - MTX IC90", fontsize=20)
g1.set_ylabel("Selection coefficient - MTX IC75", fontsize=20)
ax.set_xlim(-1,2)
ax.set_ylim(-1,2)
ax.tick_params(axis='both', labelsize=20)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
cax.tick_params(axis='y', labelsize=20)
cax.set_ylabel("Selection coefficient - MTX IC90", fontsize=20, labelpad=15)
g.axvline(0.243335428,color = 'grey', linestyle = 'dotted', linewidth = 3)
g.axvline(0.549683,color = 'black', linestyle = 'dashed', linewidth = 3)
g.axhline(0.219888764,color = 'grey', linestyle = 'dotted', linewidth = 3)
g.axhline(0.431753,color = 'black', linestyle = 'dashed', linewidth = 3)
fig.tight_layout()

fig.savefig('./Figure3A.png', dpi=300, bbox_inches='tight')
fig.show()
  
#%%
# Panel B
# Growth rate IC90 VS IC90
# Set condition
# =============================================================================
condition = 'MTX' 
# =============================================================================
# Make backup
df = data.copy()
# Drop NA
df_ok = df.dropna()

# Set data/parameters of interest depending on which condition you use
condition = "MTX"
if condition == "MTX":
    condition_num = 'MTX20'
    to_keep = ['mutant','median_MTX20', 'MTX20_max', 'MTX20_min', 'median_dgr_MTX', 'MTX_dgr_max', 'MTX_dgr_min']
    #df = data[to_keep]
    height = 1
    pal = "coolwarm"
    ex = 0.6
    IC = 'IC90'
    empty = 0.0459
    sc = 0.196
    wt = "#82A6FB"
    pp = 14
elif condition == "DMSO":
    to_keep = ['mutant','median_DMSO', 'DMSO_max', 'DMSO_min', 'median_dgr_DMSO', 'DMSO_dgr_max', 'DMSO_dgr_min']
    #df = data[to_keep]
    height = 0.2
    pal = sns.diverging_palette(135, 30, l=65, center="dark", as_cmap=True)
    ex = -0.9
    empty = 0.3075
    sc = 0.289
    wt = "#151515"
    pp = 5

#Drop control values from dataframe
df_ex = df_ok[df_ok["mutant"] != "PjDHFR"]
df_ex = df_ex[df_ex["mutant"] != "Empty"]

# Generate figure
fig, ax = plt.subplots(figsize = (6,6))

# Create scatterplot with color palette
g = sns.scatterplot(data=df_ex, x='median_MTX20', y='median_dgr_MTX', hue='median_MTX20', palette=pal, vmin=-0.5, vmax=1.5, legend=False, zorder=10, ax=ax, linewidth = 0)

# Define and plot error bars
xerror = [df_ex['MTX20_min'].to_list(),df_ex['MTX20_max'].to_list()]
yerror = [df_ex['MTX_dgr_max'].to_list(),df_ex['MTX_dgr_min'].to_list()]
g.errorbar(x=df_ex['median_MTX20'], y=df_ex['median_dgr_MTX'], xerr = xerror, yerr = yerror, color = 'black', fmt = 'o')
g.axvline

# Highlight the 'wt' data points
df_wt = df_ok[df_ok["mutant"] == "PjDHFR"]
g1 = sns.scatterplot(df_wt, x='median_MTX20', y='median_dgr_MTX', color=wt, marker='X', s=100, legend=False, zorder=10,ax=ax)

# Create regplot
g2 = sns.regplot(data=df_ex, x='median_MTX20', y='median_dgr_MTX', scatter=False, color='lightgrey',ax=ax)

# Plot outliers
out = ['G58Q', 'D110W', 'S111E', 'I196W', 'I196F']
df_out = df_ok.loc[df_ok['mutant'].isin(out)]
g3 = sns.scatterplot(data=df_out, x='median_MTX20', y='median_dgr_MTX', marker='o',ax=ax, zorder=9, s= 100, linewidth = 0, color = "green")


# Plot colorbar
norm = plt.Normalize(-0.75, 1.5)
sm = plt.cm.ScalarMappable(cmap=pal,norm = norm)
cax = fig.add_axes([ax.get_position().x1+0.1, ax.get_position().y0+0.04, 0.02, ax.get_position().height])
ax.figure.colorbar(sm, cax=cax)


# Add other plot elements
coef, p = spearmanr(df_ex['median_MTX20'], df_ex['median_dgr_MTX'])
coef = round(coef, 2)
p = round(p, 10)
g1.text(s="$\itρ$ = {}\np = {}".format(coef, p), x=ex, y=0.3, fontsize=20)
g1.set_ylabel("Growth rate - MTX IC90 (OD600$\cdot h^{-1}$)", fontsize=20)
g1.set_xlabel("Selection coefficient - MTX IC90", fontsize=20)
ax.set_xlim(-1,2)
ax.set_ylim(0,0.35)
ax.tick_params(axis='both', labelsize=20)
cax.tick_params(axis='y', labelsize=20)
cax.set_ylabel("Selection coefficient - MTX IC90", fontsize=20, labelpad=15)
g.axvline(0.243335428,color = 'grey', linestyle = 'dotted', linewidth = 2)
g.axvline(0.549683,color = 'black', linestyle = 'dashed', linewidth = 2)
g.axhline(y=empty, color="red", linestyle="-", linewidth=0.8)
g.axhline(y=sc, color="blue", linestyle="-", linewidth=0.8)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
fig.tight_layout()

fig.savefig('./Figure3B.png', dpi=300, bbox_inches='tight')
fig.show()
  
#%%
# Panel C
# Growth rate DMSO VS DMSO
# Set condition
# =============================================================================
condition = 'DMSO' 
# =============================================================================
# Make backup
df = data.copy()
# Drop NA
df_ok = df.dropna()

if condition == "MTX":
    condition_num = 'MTX20'
    to_keep = ['mutant','median_MTX20', 'MTX20_max', 'MTX20_min', 'median_dgr_MTX', 'MTX_dgr_max', 'MTX_dgr_min']
    #df = data[to_keep]
    height = 1
    pal = "coolwarm"
    ex = -0.9
    IC = 'IC90'
    empty = 0.0459
    sc = 0.196
    wt = "#82A6FB"
    pp = 14
elif condition == "DMSO":
    to_keep = ['mutant','median_DMSO', 'DMSO_max', 'DMSO_min', 'median_dgr_DMSO', 'DMSO_dgr_max', 'DMSO_dgr_min']
    #df = data[to_keep]
    height = 0.2
    pal = sns.diverging_palette(135, 30, l=65, center="dark", as_cmap=True)
    ex = -0.9
    empty = 0.3075
    sc = 0.289
    wt = "#151515"
    pp = 5


df_ex = df_ok[df_ok["mutant"] != "PjDHFR"]
df_ex = df_ex[df_ex["mutant"] != "Empty"]

fig, ax = plt.subplots(figsize = (6,6))
# Create scatterplot with the desired color palette
g = sns.scatterplot(data=df_ex, x='median_DMSO', y='median_dgr_DMSO', hue='median_DMSO', palette=pal, vmin=-0.5, vmax=1.5, legend=False, zorder=10, ax=ax, linewidth = 0)

# Highlight the 'wt' data points
df_wt = df_ok[df_ok["mutant"] == "PjDHFR"]
g1 = sns.scatterplot(df_wt, x='median_DMSO', y='median_dgr_DMSO', color=wt, marker='X', s=100, legend=False, zorder=10,ax=ax)

# Create a regression plot (regplot)
g2 = sns.regplot(data=df_ex, x='median_DMSO', y='median_dgr_DMSO', scatter=False, color='lightgrey', label='Regression Line',ax=ax)

# Create colorbar
norm = plt.Normalize(-0.08, 0.08)
sm = plt.cm.ScalarMappable(cmap=pal,norm = norm)
cax = fig.add_axes([ax.get_position().x1+0.1, ax.get_position().y0+0.04, 0.02, ax.get_position().height])
ax.figure.colorbar(sm, cax=cax)


# Add other plot elements
coef, p = spearmanr(df_ex['median_DMSO'], df_ex['median_dgr_DMSO'])
coef = round(coef, 2)
p = round(p, 2)
g1.text(s="$\itρ$ = {}\np = {}".format(coef, p), x=ex, y=0.2, fontsize=20)
g1.set_ylabel("Growth rate - DMSO (OD600$\cdot h^{-1}$)", fontsize=20)
g1.set_xlabel("Selection coefficient - DMSO", fontsize=20)
ax.set_xlim(-1,2)
ax.set_ylim(0,0.35)
ax.tick_params(axis='both', labelsize=20)
cax.tick_params(axis='y', labelsize=20)
g1.axhline(y=empty, color="red", linestyle="-", linewidth=0.8)
g1.axhline(y=sc, color="blue", linestyle="-", linewidth=0.8)
cax.set_ylabel("Selection coefficient - DMSO", fontsize=20, labelpad=15)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#g.axhline(y=empty, color="red", linestyle=":", linewidth=2)
#g.axhline(y=sc, color="blue", linestyle=":", linewidth=2)
fig.tight_layout()


fig.savefig('./Figure3C.png', dpi=300, bbox_inches='tight')
fig.show()