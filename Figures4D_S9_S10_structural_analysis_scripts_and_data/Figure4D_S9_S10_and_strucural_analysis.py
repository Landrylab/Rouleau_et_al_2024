#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 14:44:20 2023

@author: francois
"""
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
from decimal import Decimal
#%%
# Import and format data
seq ="DWQKSLTLIVALTLSRGIGLKNDLPWKLKSDMMFFSRVTSGLLVTRSTGQMNVVLMGRKTWESLPAHSRPLKNRINVVISRQEVLDLGGGAYHARSLDDALALLSQIYDSTSKIQLNRVFVIGGGELCKAAMEHSRLNRIIATVIHNEVDCDVFFPIDFRSSQSCLPWRKQDHSVLEAWVGSKVPQGKINENGFIYEFEMWIRD"
aa_list = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", 
                "A", "V", "I", "L", "M", "F", "Y", "W", "*"]
aa_list_nostop = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", 
                "A", "V", "I", "L", "M", "F", "Y", "W"]

hm_form_aa_CoS_DMSO = pd.read_csv("./total_DMSO_heatmap_format_CoS_t15_aa-prop.csv", index_col=0)
hm_form_aa_CoS_DMSO.drop('*', axis = 0, inplace = True)
hm_form_aa_CoS_DMSO_melt = pd.melt(hm_form_aa_CoS_DMSO.reset_index(), id_vars = 'index', value_vars = hm_form_aa_CoS_DMSO.columns, var_name = 'Position', value_name = 'CoS_DMSO')

hm_form_aa_CoS_MTX4 = pd.read_csv("./total_MTX4_heatmap_format_CoS_t15_aa-prop.csv", index_col=0)
hm_form_aa_CoS_MTX4.drop('*', axis = 0, inplace = True)
hm_form_aa_CoS_MTX4_melt = pd.melt(hm_form_aa_CoS_MTX4.reset_index(), id_vars = 'index', value_vars = hm_form_aa_CoS_MTX4.columns, var_name = 'Position', value_name = 'CoS_MTX4')

hm_form_aa_CoS_MTX20 = pd.read_csv("./total_MTX20_heatmap_format_CoS_t15_aa-prop.csv", index_col=0)
hm_form_aa_CoS_MTX20.drop('*', axis = 0, inplace = True)
hm_form_aa_CoS_MTX20_melt = pd.melt(hm_form_aa_CoS_MTX20.reset_index(), id_vars = 'index', value_vars = hm_form_aa_CoS_MTX20.columns, var_name = 'Position', value_name = 'CoS_MTX20')

# FOLDX DATA
pred_MTX = pd.read_csv('./Pj_DHFR_NAP_MTX_A-C_interfaces_ddgs.csv', index_col=0)
pred_MTX = pred_MTX.reindex(aa_list_nostop)
pred_MTX_melt = pd.melt(pred_MTX.reset_index(),id_vars = 'Index', value_vars = pred_MTX.columns, var_name = 'Position', value_name = 'ddg_MTX')

pred_DHF = pd.read_csv('./Pj_DHFR_NAP_FOL_A-C_interfaces_ddgs.csv', index_col=0)
pred_DHF = pred_DHF.reindex(aa_list_nostop)
pred_DHF_melt = pd.melt(pred_DHF.reset_index(), id_vars = 'Index', value_vars = pred_DHF.columns, var_name = 'Position', value_name = 'ddg_DHF')

pred_NAP = pd.read_csv('./Pj_DHFR_NAP_MTX_A-B_interfaces_ddgs.csv', index_col=0)
pred_NAP = pred_NAP.reindex(aa_list_nostop)
pred_NAP_melt = pd.melt(pred_NAP.reset_index(),id_vars = 'Index', value_vars = pred_NAP.columns, var_name = 'Position', value_name = 'ddg_NAP')

stability = pd.read_csv('./Pj_DHFR_NAP_MTX_mutations_ddgs.csv', index_col=0)
stability = stability.reindex(aa_list_nostop)
stability_melt = pd.melt(stability.reset_index(),id_vars = 'Index', value_vars = stability.columns, var_name = 'Position', value_name = 'stability')

# Now data is properly ordered
# Create complete dataframe
full_df = hm_form_aa_CoS_DMSO_melt.copy()
full_df['CoS_MTX4'] = hm_form_aa_CoS_MTX4_melt['CoS_MTX4']
full_df['CoS_MTX20'] = hm_form_aa_CoS_MTX20_melt['CoS_MTX20']
full_df['foldx_ddg_MTX'] = pred_MTX_melt['ddg_MTX']
full_df['foldx_ddg_DHF'] = pred_DHF_melt['ddg_DHF']
full_df['foldx_ddg_NAP'] = pred_NAP_melt['ddg_NAP']
full_df['stability_foldx'] = stability_melt['stability']

full_df['foldx_MTX-DHFR_ddG'] = full_df['foldx_ddg_MTX'] - full_df['foldx_ddg_DHF']

#%%
# Flexddg data for 199
flex_199_FOL = pd.read_csv("./FOL_199_results.csv", index_col=0)
flex_199_FOL = flex_199_FOL[flex_199_FOL['score_function_name'] == 'fa_talaris2014']
flex_199_FOL = flex_199_FOL[flex_199_FOL['scored_state'] == 'ddG']
flex_199_FOL = flex_199_FOL.groupby('case_name', as_index=False).median()
flex_199_FOL['aa'] = flex_199_FOL['case_name'].str.split('_').str[-1]

flex_199_MTX = pd.read_csv("./MTX_199_results.csv", index_col=0)
flex_199_MTX = flex_199_MTX[flex_199_MTX['score_function_name'] == 'fa_talaris2014']
flex_199_MTX = flex_199_MTX[flex_199_MTX['scored_state'] == 'ddG']
flex_199_MTX = flex_199_MTX.groupby('case_name', as_index=False).median()
flex_199_MTX['aa'] = flex_199_MTX['case_name'].str.split('_').str[-1]

subset = ['aa', 'total_score', 'nstruct']
flex_199_FOL = flex_199_FOL[subset]
flex_199_MTX = flex_199_MTX[subset]

flex_199_FOL = flex_199_FOL.sort_values(by='aa', key=lambda x: x.map(dict(zip(aa_list_nostop, range(len(aa_list_nostop)))))).reset_index(drop = True)
flex_199_MTX = flex_199_MTX.sort_values(by='aa', key=lambda x: x.map(dict(zip(aa_list_nostop, range(len(aa_list_nostop)))))).reset_index(drop = True)

#Stability Rosetta
flex_199_FOL_stab = pd.read_csv("./FOL_199_results.csv", index_col=0)
flex_199_FOL_stab = flex_199_FOL_stab[flex_199_FOL_stab['scored_state'] == 'mut_dG']
flex_199_FOL_stab = flex_199_FOL_stab.groupby('case_name', as_index=False).median()
flex_199_FOL_stab['aa'] = flex_199_FOL_stab['case_name'].str.split('_').str[-1]
flex_199_FOL_stab = flex_199_FOL_stab[subset]

flex_199_MTX_stab = pd.read_csv("./MTX_199_results.csv", index_col=0)
flex_199_MTX_stab = flex_199_MTX_stab[flex_199_MTX_stab['scored_state'] == 'mut_dG']
flex_199_MTX_stab = flex_199_MTX_stab.groupby('case_name', as_index=False).median()
flex_199_MTX_stab['aa'] = flex_199_MTX_stab['case_name'].str.split('_').str[-1]
flex_199_MTX_stab = flex_199_MTX_stab[subset]

full_df_199 = full_df.copy()
full_df_199 = full_df_199[full_df_199['Position'] == '199']
full_df_199.reset_index(drop = False, inplace = True)
full_df_199['rosetta_ddg_DHF'] = flex_199_FOL['total_score']
full_df_199['rosetta_ddg_MTX'] = flex_199_MTX['total_score']

#%%
# Figure S10
#Binding DHF VS MTX
xg = 'rosetta_ddg_DHF'
yg = 'rosetta_ddg_MTX'

if xg.split('_')[0] == "CoS":
    df = full_df_199.dropna()
elif yg.split('_')[0] == "CoS":
    df = full_df_199.dropna()
else:
    df = full_df_199
    
if yg == 'stability_foldx':
    full_df_199 = full_df_199[full_df['stability_foldx'] <= 7]
else:
    pass

#df = df[df['Position'] == '199']
pal = "coolwarm"
palette_cbar = ['#9ec8ff', '#dadde1', '#e0dddc', '#f9d1c1','#fbd0c0','#fccebc'
                ,'#fecebb','#ffc8b1','#ffc4aa','#ffc2a9','#ffbea2','#ffa283'
                ,'#ff9b7a','#ff9b7a','#ff8b6b','#ff8768','#ff6950','#d30023'
                ,'#d30023']

fig, ax = plt.subplots(figsize = (6,6))
r, p = spearmanr(a = df[xg], b = df[yg])
r = round(r,2)
p = round(p,4)
p = '%.1e' % Decimal(p)

plt.text(s="$\itρ$ = {}\np = {}".format(r, p), x=-0.25, y=0.4, fontsize=15)

g = sns.regplot(df, x = xg, y = yg, color = "grey", scatter = False)
g1 = sns.scatterplot(df, x = xg, y = yg, hue='CoS_MTX20', palette=palette_cbar, vmin=-0.5, vmax=1.3, legend=False, zorder=10)
norm = plt.Normalize(-0.75, 1.5)
sm = plt.cm.ScalarMappable(cmap=pal, norm = norm)
cax = fig.add_axes([ax.get_position().x1+0.1, ax.get_position().y0+0.005, 0.02, ax.get_position().height])
ax.figure.colorbar(sm, cax=cax)

#plt.title("Mutant side chain length vs Selection coefficient - IC90")
g.set_xlabel('Change in binding energy to DHF (REU)', fontsize=20)
g.set_ylabel('Change in binding energy to MTX (REU)', fontsize=20)
g.set_xlim(-0.3,0.55)
g.set_ylim(-0.3,0.55)
ax.tick_params(axis='both', labelsize=15)
cax.tick_params(axis='y', labelsize=15)
cax.set_ylabel("Selection coefficient - MTX IC90", fontsize=20, labelpad=15)
extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
g.plot([-0.3,0.55],[-0.3,0.55], linestyle = ':', color = 'black')
fig.tight_layout()
fig.savefig('./FigureS10.png', format = "png", dpi = 300, bbox_inches=extent.expanded(1.85, 1.4))


#%%
# Figure 4D
#Binding CoS_MTX vs rosetta_ddg_MTX
xg = 'CoS_MTX20'
yg = 'rosetta_ddg_MTX'
xgg = 'CoS_MTX4'

if xg.split('_')[0] == "CoS":
    df = full_df_199.dropna()
elif yg.split('_')[0] == "CoS":
    df = full_df_199.dropna()
else:
    df = full_df_199
    
if yg == 'stability_foldx':
    full_df_199 = full_df_199[full_df['stability_foldx'] <= 7]
else:
    pass
palette_cbar = ['#9ec8ff', '#dadde1', '#e0dddc', '#f9d1c1','#fbd0c0','#fccebc'
                ,'#fecebb','#ffc8b1','#ffc4aa','#ffc2a9','#ffbea2','#ffa283'
                ,'#ff9b7a','#ff9b7a','#ff8b6b','#ff8768','#ff6950','#d30023'
                ,'#d30023']
pal = "coolwarm"
fig, ax = plt.subplots(figsize = (6,6))
r20, p20 = spearmanr(a = df[xg], b = df[yg])
r20 = round(r20,2)
p20 = round(p20,3)


g = sns.regplot(df, x = xg, y = yg, color = "grey", scatter = False, ax=ax)
g1 = sns.scatterplot(df, x = xg, y = yg, hue='CoS_MTX20', palette=palette_cbar, vmin=-0.5, vmax=1.3, legend=False, ax=ax)

g.text(s="$\itρ$ = {}\np = {}".format(r20, p20), x=-0, y=0.4, fontsize=15)

norm = plt.Normalize(-0.75, 1.5)
sm = plt.cm.ScalarMappable(cmap=pal,norm = norm)
cax = fig.add_axes([ax.get_position().x1+0.1, ax.get_position().y0+0.005, 0.02, ax.get_position().height])
ax.figure.colorbar(sm, cax=cax)


# Add other plot elements
g1.set_ylabel("Change in binding energy to MTX (REU)", fontsize=20)
g1.set_xlabel("Selection coefficient - MTX IC90", fontsize=20)
ax.tick_params(axis='both', labelsize=15)
cax.tick_params(axis='y', labelsize=15)
cax.set_ylabel("Selection coefficient - MTX IC90", fontsize=20, labelpad=15)
extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
fig.tight_layout()
fig.savefig('./Figure4D.png', format = "png", dpi = 300, bbox_inches=extent.expanded(2, 1.5))
fig.show()


#%%
# Figure S9
#Binding CoS_DMSO vs stability
runcell(0)
runcell(1)
runcell(2)
xg = 'CoS_DMSO'
yg = 'stability_foldx'

if xg.split('_')[0] == "CoS":
    df = full_df.dropna()
elif yg.split('_')[0] == "CoS":
    df = full_df.dropna()
else:
    df = full_df
    
if yg == 'stability_foldx':
    df['stability_foldx'].values[df['stability_foldx'] > 5] = 5
    df['stability_foldx'].values[df['stability_foldx'] < -5] = -5
else:
    pass


plt.figure(figsize=(10,10))

r, p = spearmanr(a = df[xg], b = df[yg])
r = round(r,2)
p = round(p,len(str(p))+6)
plt.text(s="$\itρ$ = {}\np = {}".format(r, p), x=-0.2, y=-4, fontsize=15)

sns.scatterplot(df, x = xg, y = yg, color = 'black', legend=False)
sns.regplot(df, x = xg, y = yg, color = "orange", scatter = False, lowess = True)
plt.xlabel('Selection coefficient - DMSO', fontsize = 20)
plt.ylabel('Change in free energy - ΔG (kcal mol⁻¹)', fontsize = 20)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)



plt.savefig('./FigureS9.png', format = "png", dpi = 300)







