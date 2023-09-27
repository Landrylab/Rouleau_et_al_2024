#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math
import scipy.stats as stats
from scipy.stats import spearmanr
from matplotlib.ticker import FormatStrFormatter

#%%
# Validation growth curves and growth rate calculation for validations in strain ##IGA130##
# Import plate plan
source_plate = './Plate_plan_validations_IGA130.xlsx'
plate_df = pd.read_excel(source_plate, header=0, index_col=0)
plate_df.head(3)

#Import OD data
source_table = './Validations_IGA130.xlsx'
source_df = pd.read_excel(source_table, header=0, index_col=0)
source_df.rename(columns=lambda x: float(x)/3600, inplace=True)
source_df.head()

#Put data in longform
longdf = source_df.reset_index().rename(columns = {'index' : 'well'})
longdf = longdf.melt(id_vars = ['well'], var_name = 'Time (h)', value_name = 'OD')

#Annotate with layout
merged_df = longdf.merge(right=plate_df, on=['well'])
merged_df['log2_OD']= np.log2(merged_df['OD'])
merged_df = merged_df[merged_df["strain"] != "x"]
merged_df.head(5)
merged_df = merged_df[merged_df["valid_id"] != "x"]


# Define conversion function
def convert_to_uM(ug, mass_mol):
    # a simple function to convert ug/ml to uM with the molecular weigth of the compound
    vol = 0.001
    uM = round((((ug*10**-6)/(vol*mass_mol))*10**6),3)
    return uM

#%%
# Figure S5 - Plot all individual growth curves in DMSO and MTX IC90
# Define plot title order
order = pd.read_excel("./Validations_IGA130_order.xlsx")
order = order['strain']

# Create copy of imported data
merged_df_tot = merged_df
# Set values within new dataframe to ease plotting
merged_df_tot['row'] = merged_df_tot['valid_id'].str[0]
merged_df_tot['col'] = merged_df_tot['valid_id'].str[1:]
merged_df_tot = merged_df_tot.sort_values(by=["row","col"])
merged_df_tot = merged_df_tot.loc[~merged_df_tot['strain'].isin(['ScDHFR', 'PjDHFR', 'Empty', 'pKB44'])]
row_ord = ["A","B","C","D","E","F","G","H"]
col_ord = ["1","2","3","4","5","6"]

# Convert g/mL units into uM
# MTX molecular mass: 454.44g/mol
merged_df_tot['concentration_uM'] = merged_df_tot[merged_df_tot.drug=='MTX']['concentration'].apply(lambda x: convert_to_uM(x, 454.44))
merged_df_tot = merged_df_tot.fillna(0)

# Plot data
grid = sns.FacetGrid(data=merged_df_tot[merged_df_tot['Time (h)'] <=40],
                     col="col",
                     row = "row",
                     hue="concentration_uM", palette='rocket')

# Other graph information
grid.map(sns.lineplot, 'Time (h)', 'OD')
grid.set_axis_labels('Time (h)', '$\mathregular{OD_{600}}$')
grid.add_legend(title='MTX concentration (μM)')
grid.fig.subplots_adjust(top=0.9)
grid.tight_layout()
titles = order
for ax, title in zip(grid.axes.flatten(), titles):
    ax.set_title(title)
grid.tight_layout()

# Save figure
plt.savefig('./FigureS6.png', format='png',dpi=400)

#%%
# Get growth rate from derivative

# Generate intermediate dataframe for data storage
temp_df = pd.DataFrame()
dg_df = pd.DataFrame(columns = ["well", "mutant", "dgr"])

# =============================================================================
concentration = 0 # 0 for DMSO, 20 for IC90
# =============================================================================

# Compute data to get growth rate and remove unecessary values
start_df = merged_df[merged_df["concentration"] == concentration]
start_df = start_df[start_df["valid_id"] != "pKB44"]

# Remove the first hour, caused issues between machines
crop = ['0','0.25','0.5','0.75']
start_df = start_df[~start_df['Time (h)'].isin(crop)]

# Get the information for all wells in plate plan
wells = start_df.well
wells = list(dict.fromkeys(wells))
index = 0
for well in wells:
    temp_df = start_df[start_df["well"] == well]
    mutant = temp_df.iloc[0,3]
    temp_df["dg"] = temp_df["OD"].pct_change()
    top5 = temp_df["dg"].nlargest(5)
    top5_med = (top5.median())*4
    list_val = [well, mutant, top5_med]
    dg_df.loc[index] = list_val
    temp_df = pd.DataFrame()
    index = index + 1

# Output maximum growth rate median for each validation mutant in a given condition
dg_df_group = dg_df.groupby(by='mutant', as_index = False).median()
#%%
# Get growth rates for all individual replicates to compute error bars in other figures
print_out = list(dict.fromkeys(dg_df.mutant))
for mutants in print_out:
    dg_df_print = dg_df[dg_df['mutant'] == mutants]
    mutant_array = dg_df_print['dgr'].to_list()
    print(mutants)
    print(str(mutant_array[0])+'\t'+str(mutant_array[1])+'\t'+str(mutant_array[2]))

#%%
# Validation growth curves and growth rate calculation for validations in strain ##FDR0001##
runcell(0)
## Growth curves were done with and without estradiol, specify conditions here ##
estradiol = 100 # 0 for no estradiol, 100 for with estradiol
concentration = 20 # DMSO for no drug, 20 for IC90

# Import plate plan
source_plate = './Plate_plan_validations_FDR0001.xlsx'
plate_df = pd.read_excel(source_plate, header=0, index_col=0)
plate_df.head(3)

#Import OD data
source_table = './Validations_FDR0001.xlsx'
source_df = pd.read_excel(source_table, header=0, index_col=0)
source_df.rename(columns=lambda x: float(x)/3600, inplace=True)
source_df.head()

#Put data in longform
longdf = source_df.reset_index().rename(columns = {'index' : 'well'})
longdf = longdf.melt(id_vars = ['well'], var_name = 'Time (h)', value_name = 'OD')

#Annotate with layout
merged_df = longdf.merge(right=plate_df, on=['well'])
merged_df['log2_OD']= np.log2(merged_df['OD'])
merged_df = merged_df[merged_df["strain"] != "x"]
merged_df.head(5)
merged_df = merged_df[merged_df["valid_id"] != "x"]


# Define conversion function
def convert_to_uM(ug, mass_mol):
    # a simple function to convert ug/ml to uM with the molecular weigth of the compound
    vol = 0.001
    uM = round((((ug*10**-6)/(vol*mass_mol))*10**6),3)
    return uM

# Generate temporary variables
temp_df = pd.DataFrame()
dg_df = pd.DataFrame(columns = ["well", "mutant", "dgr"])
index = 0

## Set desired parameters
start_df = merged_df[merged_df["concentration"] == concentration]
start_df = start_df[start_df["estradiol"] == estradiol]

## Crop initial values (machine used sometimes missread the first hour)
crop = ['0','0.25','0.5','0.75']
start_df = start_df[~start_df['Time (h)'].isin(crop)]


# Get growth rates
wells = start_df.well
wells = list(dict.fromkeys(wells))
for well in wells:
    temp_df = start_df[start_df["well"] == well]
    mutant = temp_df.iloc[0,3]
    temp_df["dg"] = temp_df["OD"].pct_change()
    top5 = temp_df["dg"].nlargest(5)
    top5_med = (top5.median())*4
    list_val = [well, mutant, top5_med]
    dg_df.loc[index] = list_val
    temp_df = pd.DataFrame()
    index = index + 1

# Get median
dg_df_group = dg_df.groupby(by='mutant', as_index = False).median()

# Get each value independantly
## Print values to transfer them to error bar file
print_out = list(dict.fromkeys(dg_df.mutant))
for mutants in print_out:
    dg_df_print = dg_df[dg_df['mutant'] == mutants]
    mutant_array = dg_df_print['dgr'].to_list()
    print(mutants)
    print(str(mutant_array[0])+'\t'+str(mutant_array[1])+'\t'+str(mutant_array[2])+'\t'+str(mutant_array[3]))

#%%
# Load information outputed from previous cell and arranged in table manually
if estradiol == 100:
    back = pd.read_excel("./CoS_dgr_replicates_error_FDR_full_estra.xlsx", index_col = 0)
elif estradiol == 0:
    back = pd.read_excel("./CoS_dgr_replicates_error_FDR_full.xlsx", index_col = 0)
data_error = back.copy()
data_error[['MTX_dgr_max','MTX_dgr_min','DMSO_dgr_max','DMSO_dgr_min']] = np.nan
data_error = data_error.reset_index(drop = True)


#%%
# Calculate data for error bars
medians = ['median_dgr-MTX', 'median_dgr-DMSO']
for row in data_error.index:
    data_mutant = data_error.iloc[row,]
    MTXdgr_median = data_mutant.at['median_dgr-MTX']
    array_dgr_MTX = [data_mutant.at['dgr1_MTX'],data_mutant.at['dgr2_MTX'],data_mutant.at['dgr3_MTX'],data_mutant.at['dgr4_MTX']]
    
    DMSOdgr_median = data_mutant.at['median_dgr-DMSO']
    array_dgr_DMSO = [data_mutant.at['dgr1_DMSO'],data_mutant.at['dgr2_DMSO'],data_mutant.at['dgr3_DMSO'],data_mutant.at['dgr4_DMSO']]

    array_array = [array_dgr_MTX,array_dgr_DMSO]
    
    for array in array_array:
        if array == array_dgr_MTX:
            median = MTXdgr_median
            for value in array:
                if value > median:
                    up = value - median
                    data_error.at[row,'MTX_dgr_max'] = up
                if value < median:
                    down = median - value
                    data_error.at[row,'MTX_dgr_min'] = down
        if array == array_dgr_DMSO:
            median = DMSOdgr_median
            for value in array:
                if value > median:
                    up = value - median
                    data_error.at[row,'DMSO_dgr_max'] = up
                if value < median:
                    down = median - value
                    data_error.at[row,'DMSO_dgr_min'] = down
  
  
data_error = data_error.fillna(0)
#%%
# Construct plot with error bars
## Load data
source_plate = './FR001_growthrate_vs_DMS.xlsx'
df = pd.read_excel(source_plate, header=0)
graph = pd.DataFrame()

# Drop NA & subset data depending on estradiol concentration
df_ok = df.copy()
if estradiol == 100:
    subset = df_ok[df_ok["estradiol"] == 100]
elif estradiol == 0:
    subset = df_ok[df_ok["estradiol"] == 0]

# Rearrange data for plot
complementation = subset[subset["drug"] == "DMSO"]
complementation.reset_index(inplace = True, drop = True)
resistance = subset[subset["drug"] == "MTX"]
resistance.reset_index(inplace = True, drop = True)
graph["mutant"] = complementation["mutant"]
graph["complementation"] = complementation["dgr"]
graph["resistance"] = resistance["dgr"]
graph["CoS"] = resistance["CoS20"]
ctrls = ["Empty", "PjDHFR", 'ScDHFR', "mDHFR"]
df_wt = graph.loc[graph['mutant'].isin(ctrls)]
graph = graph.loc[~graph['mutant'].isin(ctrls)]

# Construct plot
fig, ax = plt.subplots(figsize = (8,8))
g1 = ax.scatter(data = graph, x = "complementation", y = "resistance", 
                    c=graph["CoS"], cmap="coolwarm", vmin=-0.75, vmax=1.5 
                    ,zorder=10,linewidth = 0)


xerror = [data_error['DMSO_dgr_min'].to_list(),data_error['DMSO_dgr_min'].to_list()]
yerror = [data_error['MTX_dgr_min'].to_list(),data_error['MTX_dgr_min'].to_list()]
ax.errorbar(x=data_error['median_dgr-DMSO'], y=data_error['median_dgr-MTX'], xerr = xerror, yerr = yerror, color = 'black', fmt = 'o')


g2 = sns.scatterplot(df_wt, x='complementation', y='resistance', color=['#d55e00',"#0072b2","#f0e442","#009e73"],edgecolor="black", marker='X', s=100, legend=False, zorder=10,ax=ax)

norm = plt.Normalize(-0.75, 1.5)
sm = plt.cm.ScalarMappable(cmap="coolwarm",norm = norm)
cax = fig.add_axes([ax.get_position().x1+0.1, ax.get_position().y0+0.04, 0.02, ax.get_position().height])
ax.figure.colorbar(sm, cax=cax)

# Add other plot elements
coef, p = spearmanr(graph['complementation'], graph['resistance'])
coef = round(coef, 2)
if estradiol == 100:
    p = round(p, 8)
elif estradiol == 0:
    p = round(p, 3)
ax.text(s="ρ = {}\np = {}".format(coef, p), x=0.02, y=0.25, fontsize=20)
ax.set_ylabel("Growth rate - MTX IC90 (OD600$\cdot h^{-1}$)", fontsize=20)
ax.set_xlabel("Growth rate - DMSO (OD600$\cdot h^{-1}$)", fontsize=20)
ax.set_xlim(0,0.35)
ax.set_ylim(0,0.35)
ax.tick_params(axis='both', labelsize=20)
cax.tick_params(axis='y', labelsize=20)
ax.axline([ax.get_xlim()[0], ax.get_ylim()[0]], [ax.get_xlim()[1], ax.get_ylim()[1]], color = "black", linestyle = ":")
fig.tight_layout()
if estradiol == 100:
    plt.savefig("./FigureS7B.png", format = "png", dpi = 300, bbox_inches='tight')
elif estradiol == 0:
    plt.savefig("./FigureS7A.png", format = "png", dpi = 300, bbox_inches='tight')