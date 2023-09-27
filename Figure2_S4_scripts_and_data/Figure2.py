#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 14:44:20 2023

@author: francois
"""
# Figure 2. Pannels are to be assembled individually.
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
# Run data aggregation for all conditions
conditions = ['DMSO','MTX4','MTX20']
for condition in conditions:
    longform = pd.read_csv("../Longform_dataframes/total_"+condition+"_longform_CoS_t15.csv")
    aa_list = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", 
                    "A", "V", "I", "L", "M", "F", "Y", "W", "*"]
    seq ="DWQKSLTLIVALTLSRGIGLKNDLPWKLKSDMMFFSRVTSGLLVTRSTGQMNVVLMGRKTWESLPAHSRPLKNRINVVISRQEVLDLGGGAYHARSLDDALALLSQIYDSTSKIQLNRVFVIGGGELCKAAMEHSRLNRIIATVIHNEVDCDVFFPIDFRSSQSCLPWRKQDHSVLEAWVGSKVPQGKINENGFIYEFEMWIRD"
    positions = list(range(2,206))
    
    # Run cell to load condition data
    
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
    
    
    #t,p=stats.ttest_ind(longform_silent_melt[condition+'_CoS'], longform['DMSO-1_CoS'], equal_var = False)
    silent_array = longform_silent_melt[condition+'_CoS']
    
    longform_aa = pd.DataFrame(columns=['aa_WT','position','aa_mut','p-values','Welch_t-test','codons_num','median_CoS', 'mut_code', 'nature'])
    index = 0
    
    for posit in positions:
        mutant_df_master = longform[longform['position'] == posit]
        for aa in aa_list:
            aa_df = mutant_df_master[mutant_df_master['aa_mut'] == aa]
            if len(aa_df) > 0:
                nature = str(aa_df['nature'].iloc[0])
                aa_wt = seq[posit-2]
                mutant_df_array = aa_df[cols_to_keep2]
                mutant_df_median = aa_df['median_CoS'].median()
                mutant_array = mutant_df_array.melt(value_vars = cols_to_keep2, ignore_index = True)
                mutant_array = mutant_array['value']
                t, p = stats.ttest_ind(silent_array, mutant_array, equal_var = False, alternative='less')
                longform_aa.at[index,'aa_WT'] = aa_wt
                longform_aa.at[index,'position'] = posit
                longform_aa.at[index,'aa_mut'] = aa
                longform_aa.at[index,'p-values'] = p
                longform_aa.at[index,'Welch_t-test'] = t
                longform_aa.at[index,'codons_num'] = len(mutant_array)/3
                longform_aa.at[index,'median_CoS'] = mutant_df_median
                longform_aa.at[index,'mut_code'] = aa_wt + str(posit) + aa
                longform_aa.at[index,'nature'] = nature
                index = index + 1
            else:
                pass
    longform_aa = longform_aa.astype({'aa_WT': str, 
                                      'position': 'int64',
                                      'aa_mut': 'object',
                                      'p-values': 'float64',
                                      'Welch_t-test': 'float64',
                                      'codons_num': 'float64',
                                      'median_CoS': 'float64',
                                      'nature': str}, copy = True)
    if condition == "DMSO":
        longform_aa_DMSO = longform_aa.copy()
    elif condition == "MTX4":
        longform_aa_MTX4 = longform_aa.copy()
    elif condition == "MTX20":
        longform_aa_MTX20 = longform_aa.copy()

#%%
# Panel A
# Set custom parameters for plotting
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

# Load data for relevant conditions
DMSO_dat = longform_aa_DMSO.copy()
MTX20_dat = longform_aa_MTX4.copy()

# Generate dataframe for plot
paired = pd.DataFrame()
cols_sub = ["aa_WT", "position","aa_mut"]
paired[cols_sub] = DMSO_dat[cols_sub]
paired["DMSO"] = DMSO_dat["median_CoS"]
paired["MTX"] = MTX20_dat["median_CoS"]
paired["Nature"] = MTX20_dat["nature"]
paired = paired[paired['Nature']!='WT']
mapping = {'silent': 'Synonymous', 'substitution': 'Non-synonymous', 'stop': 'Nonsense'}
# Replace the values in the 'Nature' column and update the column name
paired['Nature'] = paired['Nature'].replace(mapping)
paired = paired.rename(columns={'Nature': 'Type of mutation'})

# Calculate spearman rho
rho, p_value = stats.spearmanr(paired["DMSO"], paired["MTX"])
rho = round(rho,2)
p_value = round(p_value,12)

# Plot information
hue_order = ['Non-synonymous','Synonymous','Nonsense']
colors = {'Synonymous': '#f0e442', 'Nonsense': '#d55e00', 'Non-synonymous': '#0072b2'}

# Plot
g = sns.jointplot(data=paired.sort_values('Type of mutation', key=np.vectorize(hue_order.index))
              , x='DMSO', y='MTX', kind='scatter', hue = 'Type of mutation', hue_order=hue_order, palette=colors, 
              height = 8, linewidth = 0.1)
g.ax_marg_x.remove()
g.ax_marg_y.remove()
# Set x-axis and y-axis limits and titles
plt.xlim(-1, 2)
plt.ylim(-1, 2)
plt.text(1,1.2,'$\itρ$ = {} \np = {}'.format(rho, p_value), fontsize=18)
plt.legend(framealpha = 0, fontsize = 15, borderaxespad = -0.5,
          markerscale = 1.5, handletextpad = -0.2, loc = (0.53,0))
plt.xlabel('Selection coefficient - DMSO', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.axhline(0.219888764,color = 'grey', linestyle = 'dotted', linewidth = 2)
plt.axhline(0.431753,color = 'black', linestyle = 'dashed', linewidth = 2)
plt.ylabel('Selection coefficient - MTX IC75', fontsize=20)
plt.tight_layout()
plt.savefig('./Figure2A_font.png',format = 'png', dpi=300)
plt.show()
#%%
# Panel B

# Set custom parameters for plotting
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

# Load data for relevant conditions
DMSO_dat = longform_aa_DMSO.copy()
MTX20_dat = longform_aa_MTX20.copy()

# Generate dataframe for plot
paired = pd.DataFrame()
cols_sub = ["aa_WT", "position","aa_mut"]
paired[cols_sub] = DMSO_dat[cols_sub]
paired["DMSO"] = DMSO_dat["median_CoS"]
paired["MTX"] = MTX20_dat["median_CoS"]
paired["Nature"] = MTX20_dat["nature"]
paired = paired[paired['Nature']!='WT']
mapping = {'silent': 'Synonymous', 'substitution': 'Non-synonymous', 'stop': 'Nonsense'}
# Replace the values in the 'Nature' column and update the column name
paired['Nature'] = paired['Nature'].replace(mapping)
paired = paired.rename(columns={'Nature': 'Type of mutation'})

# Calculate spearman rho
rho, p_value = stats.spearmanr(paired["DMSO"], paired["MTX"])
rho = round(rho,2)
p_value = round(p_value,20)

# Plot information
hue_order = ['Non-synonymous','Synonymous','Nonsense']
colors = {'Synonymous': '#f0e442', 'Nonsense': '#d55e00', 'Non-synonymous': '#0072b2'}

# Plot
g = sns.jointplot(data=paired.sort_values('Type of mutation', key=np.vectorize(hue_order.index))
              , x='DMSO', y='MTX', kind='scatter', hue = 'Type of mutation', hue_order=hue_order, palette=colors,
              height = 8, linewidth = 0.1)
g.ax_marg_x.remove()
g.ax_marg_y.remove()
# Set x-axis and y-axis limits and titles
plt.xlim(-1, 2)
plt.ylim(-1, 2)
plt.text(1,1.2,'$\itρ$ = {} \np = {}'.format(rho, p_value), fontsize=18)
plt.legend(framealpha = 0, fontsize = 15, borderaxespad = -0.5,
          markerscale = 1.5, handletextpad = -0.2, loc = (0.53,0))
plt.xlabel('Selection coefficient - DMSO', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.axhline(0.549683,color = 'black', linestyle = 'dashed', linewidth = 2)
plt.axhline(0.243335428,color = 'grey', linestyle = 'dotted', linewidth = 2)
plt.ylabel('Selection coefficient - MTX IC90', fontsize=20)
plt.tight_layout()
plt.savefig('./Figure2B_font.png',format = 'png', dpi=300)
plt.show()
#%%
# Panel D

# Set custom parameters for plotting
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

# Load data for relevant conditions
DMSO_dat = longform_aa_MTX4.copy()
MTX20_dat = longform_aa_MTX20.copy()

# Generate dataframe for plot
paired = pd.DataFrame()
cols_sub = ["aa_WT", "position","aa_mut"]
paired[cols_sub] = DMSO_dat[cols_sub]
paired["DMSO"] = DMSO_dat["median_CoS"]
paired["MTX"] = MTX20_dat["median_CoS"]
paired["Nature"] = MTX20_dat["nature"]
paired = paired[paired['Nature']!='WT']
mapping = {'silent': 'Synonymous', 'substitution': 'Non-synonymous', 'stop': 'Nonsense'}
# Replace the values in the 'Nature' column and update the column name
paired['Nature'] = paired['Nature'].replace(mapping)
paired = paired.rename(columns={'Nature': 'Type of mutation'})

# Calculate spearman rho
rho, p_value = stats.spearmanr(paired["DMSO"], paired["MTX"])
rho = round(rho,2)
p_value = round(p_value,50)

# Plot information
hue_order = ['Non-synonymous','Synonymous','Nonsense']
colors = {'Synonymous': '#f0e442', 'Nonsense': '#d55e00', 'Non-synonymous': '#0072b2'}

# Plot
g = sns.jointplot(data=paired.sort_values('Type of mutation', key=np.vectorize(hue_order.index))
              , x='MTX', y='DMSO', kind='scatter', hue = 'Type of mutation', hue_order=hue_order, palette=colors
              , height = 8, linewidth = 0.1)
g.ax_marg_x.remove()
g.ax_marg_y.remove()
# Set x-axis and y-axis limits and titles
plt.xlim(-1, 2)
plt.ylim(-1, 2)
plt.text(-0.8,1.5,'$\itρ$ = {} \np = {}'.format(rho, p_value), fontsize=18)
plt.legend(framealpha = 0, fontsize = 15, borderaxespad = -0.5,
          markerscale = 1.5, handletextpad = -0.2, loc = (0.53,0))
plt.xlabel('Selection coefficient - MTX IC90', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.axvline(0.243335428,color = 'grey', linestyle = 'dotted', linewidth = 2)
plt.axvline(0.549683,color = 'black', linestyle = 'dashed', linewidth = 2)
plt.axhline(0.219888764,color = 'grey', linestyle = 'dotted', linewidth = 2)
plt.axhline(0.431753,color = 'black', linestyle = 'dashed', linewidth = 2)
plt.ylabel('Selection coefficient - MTX IC75', fontsize=20)
plt.tight_layout()
plt.savefig('./Figure2C_font.png',format = 'png', dpi=300)
plt.show()

#%%
runcell(0)

seq ="DWQKSLTLIVALTLSRGIGLKNDLPWKLKSDMMFFSRVTSGLLVTRSTGQMNVVLMGRKTWESLPAHSRPLKNRINVVISRQEVLDLGGGAYHARSLDDALALLSQIYDSTSKIQLNRVFVIGGGELCKAAMEHSRLNRIIATVIHNEVDCDVFFPIDFRSSQSCLPWRKQDHSVLEAWVGSKVPQGKINENGFIYEFEMWIRD"
aa_list = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", 
                "A", "V", "I", "L", "M", "F", "Y", "W", "*"]


hm_form_aa_CoS_MTX20 = pd.read_csv("./total_MTX20_heatmap_format_CoS_t15_aa-prop.csv", index_col=0)

fig, ax = plt.subplots(2, 2, figsize=(24, 9), layout = 'tight',sharex = False ,gridspec_kw={'width_ratios': [25,0.5],'height_ratios': [1,0.25]})


g3 = sns.heatmap(hm_form_aa_CoS_MTX20, square=False, cmap='coolwarm', cbar_ax=ax[0,1], ax=ax[0,0]
    ,vmin = -0.75, vmax = 1.5
    )

ax[0,0].tick_params(axis='x', which='major', bottom=False, labelbottom=False)
ax[0,0].tick_params(axis='y', which='major', labelsize=20, rotation = 0.25)
ax[0,0].set_ylabel("Amino acid", fontsize=25)
#ax[0,0].set_xlabel("Position", fontsize=25)

ax[0,1].tick_params(labelsize=20)
ax[0,1].set_ylabel("Selection coefficient - MTX IC90", fontsize=25, labelpad=15)

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
, robust = True, linecolor = 'black', linewidths = 1, ax=ax[1,0])

ax[1,0].set_xlabel("Position", fontsize=25)
ax[1,0].set_ylabel("Contact", fontsize=25)
ax[1,0].tick_params(axis='y', which='major', labelsize=20, labelrotation = 30)
ax[1,0].tick_params(axis='x', which='major', labelsize=18)
ax[1,1].set_visible(False)



plt.tight_layout()
plt.savefig('./Figure2D_E.png', dpi=200, bbox_inches="tight", transparent = False)