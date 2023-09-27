#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 14:44:20 2023

@author: francois
"""
# Supp heatmap DMSO
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
#%%
runcell(0)

seq ="DWQKSLTLIVALTLSRGIGLKNDLPWKLKSDMMFFSRVTSGLLVTRSTGQMNVVLMGRKTWESLPAHSRPLKNRINVVISRQEVLDLGGGAYHARSLDDALALLSQIYDSTSKIQLNRVFVIGGGELCKAAMEHSRLNRIIATVIHNEVDCDVFFPIDFRSSQSCLPWRKQDHSVLEAWVGSKVPQGKINENGFIYEFEMWIRD"
aa_list = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", 
                "A", "V", "I", "L", "M", "F", "Y", "W", "*"]

cust = sns.diverging_palette(135, 30, l=65, center="dark", as_cmap=True)

hm_form_aa_CoS_DMSO = pd.read_csv("./total_DMSO_heatmap_format_CoS_t15_aa-prop.csv", index_col=0)
hm_form_aa_CoS_MTX4 = pd.read_csv("./total_MTX4_heatmap_format_CoS_t15_aa-prop.csv", index_col=0)

fig, ax = plt.subplots(3, 2, figsize=(24, 16.2), layout = 'tight',sharex = False ,gridspec_kw={'width_ratios': [25,0.5],'height_ratios': [1,1,0.25]})


g1 = sns.heatmap(hm_form_aa_CoS_DMSO, square=False, cmap=cust, cbar_ax=ax[0,1], ax=ax[0,0]
    ,vmin = -0.08, vmax = 0.08
    )

g2 = sns.heatmap(hm_form_aa_CoS_MTX4, square=False, cmap='coolwarm', cbar_ax=ax[1,1], ax=ax[1,0]
    ,vmin = -0.75, vmax = 1.5
    )


ax[0,0].tick_params(axis='x', which='major', bottom=False, labelbottom=False)
ax[0,0].tick_params(axis='y', which='major', labelsize=20, rotation = 0.25)
ax[0,0].set_ylabel("Amino acid", fontsize=25)


ax[1,0].tick_params(axis='x', which='major', bottom=False, labelbottom=False)
ax[1,0].tick_params(axis='y', which='major', labelsize=20, rotation = 0.25)
ax[1,0].set_ylabel("Amino acid", fontsize=25)

ax[0,1].tick_params(labelsize=20)
ax[0,1].set_ylabel("Selection coefficient - DMSO", fontsize=25, labelpad=15)
ax[1,1].tick_params(labelsize=20)
ax[1,1].set_ylabel("Selection coefficient - IC75", fontsize=25, labelpad=15)

for x in range(0,len(seq)):
        x_pos = x+0.5
        # set x coordinates of the annotation. 0.5 increments place the dots in the middle of the squares
        codon_pos = list()
        codon_n = codon_pos.append(x)
        # get the codon number in integer format
        seq = list(seq)
        # get the dna sequence of this codon number
        y_pos = float(aa_list.index(seq[x]))+0.5
        #print(y_pos, seq[x])
        # get y coordinates based on the order in in which amino acids are (the df index or a custom order if it has been changed)
        # 0.5 increments again to be in the center of the squares
        ax[0,0].plot(x_pos, y_pos, 'ro', markersize = 3.5)
        ax[1,0].plot(x_pos, y_pos, 'ro', markersize = 3.5)
        
pos_prop = pd.read_excel("./Pos_properties.xlsx")
#crest_r
dist = pd.melt(pos_prop, value_vars=['position'], id_vars=['min_dist_fol','min_dist_mtx','min_dist_NADP'
    #                                                       ,'allostery'
                                                           ])
dist_t = dist.transpose()
column_names_row = dist_t.iloc[4]
dist_t = dist_t.set_axis(column_names_row, axis=1)
rows_to_d = ['variable','value']
r_names_row = ['Distance to folate','Distance to MTX','Distance to NADPH'
        #       ,'Allostery'
               ]
dist_t = dist_t.drop(rows_to_d)
dist_t = dist_t.set_axis(r_names_row, axis=0)
dist_t = dist_t.astype(float)

df_new = dist_t.copy()

# Generate 'Contact to MTX' column
df_new.loc['DHF'] = (dist_t.loc['Distance to folate'] <= 8).astype(int)

# Generate 'Contact to MTX' column
df_new.loc['MTX'] = (dist_t.loc['Distance to MTX'] <= 8).astype(int)

# Generate 'Contact to NADPH' column
df_new.loc['NADPH'] = (dist_t.loc['Distance to NADPH'] <= 8).astype(int)

rows_to_d2 = ['Distance to folate','Distance to MTX','Distance to NADPH']
df_new = df_new.drop(rows_to_d2)

df_new2 = df_new.copy()
df_new2.loc['NADPH', df_new2.loc['NADPH'] == 1] = 2
df_new2.loc['DHF', df_new2.loc['DHF'] == 1] = 3
desired_order = ['DHF','MTX','NADPH']

# Reorder the columns based on the desired_order list
df_ordered = df_new2.reindex(index=desired_order)
df_ordered = df_ordered.drop(columns=[df_ordered.columns[0], df_ordered.columns[-1]])

allblack = ['#FFFFFF', '#151515', '#151515', '#151515']
normal = ['#FFFFFF', '#ffff00', '#1a8d1a', '#0000FF']

g4 = sns.heatmap(df_ordered, square=False, cmap=normal, cbar=False
, robust = True, linecolor = 'black', linewidths = 1, ax=ax[2,0])

ax[2,0].set_xlabel("Position", fontsize=25)
ax[2,0].set_ylabel("Contact", fontsize=25)
ax[2,0].tick_params(axis='y', which='major', labelsize=20, labelrotation = 30)
ax[2,0].tick_params(axis='x', which='major', labelsize=18)
ax[2,1].set_visible(False)



plt.tight_layout()
plt.savefig('./FigureS4AB.png', dpi=200, bbox_inches="tight", transparent = False)