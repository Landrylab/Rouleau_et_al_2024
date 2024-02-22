#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 08:46:31 2024

@author: francois
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import spearmanr
#%%
# Load raw data from cytometry assay
source_plate = './PjDHFR_Mutant_raw_cyto_mEGFP.csv'
plate_df = pd.read_csv(source_plate, header=0, index_col=0)
df = plate_df.copy()
# Drop irrelavant information
df = df[['id', 'rep', 'FSC.HLin', 'GRN.B.HLin']]
# Correct linear fluorescence intensity based on cell size and transform into log10
df["log_greenlin_corr"] = np.log(df['GRN.B.HLin']/df['FSC.HLin'])
# Get wild-type median
df_wt = df[df['id'] == "PjDHFR-mEGFP"]
median_wt = df_wt["log_greenlin_corr"].median()
#%%
# Measure linear increase in fluorescence
test = df.copy()
test = test.groupby(by = 'id', as_index = False).mean()
test["inc"] = (-3.40662/test["log_greenlin_corr"])
#%%
# Generate figures
# Boxplot with all mutants
order = ["mEGFP", "PjDHFR-mEGFP", "B1", "C1", "D2", "A3", "D3", "C4", "D4", "E5", "E6", "D6", "Empty"]

fig, ax = plt.subplots(figsize = (6,6))
g1 = sns.boxplot(data = df, x = 'id', y = "log_greenlin_corr", order = order, color = 'grey')
labels = ["mEGFP", "PjDHFR\nmEGFP", "I10H_res", "I10R_sens", 
          "M33P_res", "G58Q_sens", "D110W_sens", "G124S_res", 
          "G124L_sens", "I142G_sens", "F199E_res", "F199P_res", "Empty"]
plt.axhline(median_wt, color = 'black', linewidth = 2.5, linestyle = ":")
ax.set_xticklabels(labels)
ax.set_ylabel("Log relative fluorescence (AU) \n Corrected for cell size", fontsize=12)
ax.set_xlabel("Mutant", fontsize=12)
ax.tick_params(axis='both', labelsize=12)
plt.xticks(rotation = 90)

plt.savefig("./FigureS12A.png", dpi = 200, format = "png", bbox_inches='tight')
#%%
# Regplot with median and growth rate in MTX
source_plate1 = './Sel-coef_vs_Cyto_vs_dgr.xlsx'
plate_df1 = pd.read_excel(source_plate1, header=0, index_col=0)
df2 = plate_df1.copy()
df2 = df2[df2["Status"] != 'none']

fig, ax = plt.subplots(figsize = (6,6))
g1 = sns.regplot(data = df2, x = 'dgr_MTX', y = "median_fluo_corr", scatter=False, color = "grey")
g2 = sns.scatterplot(data = df2, x = 'dgr_MTX', y = "median_fluo_corr", hue = "Status")
ax.set_ylabel("Log relative fluorescence (AU) \n Corrected for cell size", fontsize=12)
ax.set_xlabel("Growth rate - MTX IC90 (OD600$\cdot h^{-1}$)", fontsize=12)
ax.tick_params(axis='both', labelsize=12)
leg = ax.legend(loc = "upper left", title = "Status", fontsize = 12)
leg.set_title('Status',prop={'size':12})
coef, p = spearmanr(df2['dgr_MTX'], df2['median_fluo_corr'])
coef = round(coef, 2)
p = round(p, 2)
ax.text(s="ρ = {}\np = {}".format(coef, p), x=0.035, y=-4.75, fontsize=12)

plt.savefig("./FigureS12B.png", dpi = 200, format = "png", bbox_inches='tight')











