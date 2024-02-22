#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import spearmanr

#%%
# Validation growth curves and growth rate calculation for validations in strain IGA130 or FDR0001
strain = "FDR0001"

# Import plate plan
source_plate = 'Plate_plan_'+strain+'_valids_total.xlsx'
plate_df = pd.read_excel(source_plate, header=0, index_col=0)
plate_df.head(3)

#Import OD data
source_table = './Validations_'+strain+'_total.xlsx'
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
order = pd.read_excel(".Validations_"+strain+"_order.xlsx")
order = order['strain']

# Create copy of imported data
merged_df_tot = merged_df
# Set values within new dataframe to ease plotting
merged_df_tot['row'] = merged_df_tot['valid_id'].str[0]
merged_df_tot['col'] = merged_df_tot['valid_id'].str[1:]
merged_df_tot = merged_df_tot.sort_values(by=["row","col"])
merged_df_tot = merged_df_tot.loc[~merged_df_tot['strain'].isin(['ScDHFR', 'PjDHFR', 'Empty', 'mDHFR'])]
row_ord = ["A","B","C","D","E","F","G","H"]
col_ord = ["1","2","3","4","5","6"]

# Convert g/mL units into uM
# MTX molecular mass: 454.44g/mol
merged_df_tot['concentration_uM'] = merged_df_tot[merged_df_tot.drug=='MTX']['concentration'].apply(lambda x: convert_to_uM(x, 454.44))
merged_df_tot = merged_df_tot.fillna(0)
merged_df_tot = merged_df_tot[merged_df_tot["estradiol"] == 0]

# Plot data
grid = sns.FacetGrid(data=merged_df_tot[merged_df_tot['Time (h)'] <=40],
                     col="col",
                     row = "row",
                     hue="drug", palette='rocket')

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
if strain == "IGA130":
    plt.savefig('./FigureS7.png', format='png',dpi=400)
else:
    plt.show()

#%%
# Get growth rate from derivative

# Generate intermediate dataframe for data storage
temp_df = pd.DataFrame()
dg_df = pd.DataFrame(columns = ["well", "mutant", "dgr"])

# =============================================================================
concentration = 20 # 0 for DMSO, 20 for IC90
estra = 0 # 0 for no estra, 100 for estra
# =============================================================================

# Compute data to get growth rate and remove unecessary values
start_df = merged_df[merged_df["concentration"] == concentration]
start_df = start_df[start_df["estradiol"] == estra]
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
## Growth curves were done with and without estradiol, specify conditions here ##
estradiol = 100 # 0 for no estradiol, 100 for with estradiol
concentration = 20 # DMSO for no drug, 20 for IC90

# Import plate plan
source_plate = './Plate_plan_'+strain+'_valids_total.xlsx'
plate_df = pd.read_excel(source_plate, header=0, index_col=0)
plate_df.head(3)

#Import OD data
source_table = './Validations_'+strain+'_total.xlsx'
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
    if len(mutant_array) > 3:
        print(str(mutant_array[0])+'\t'+str(mutant_array[1])+'\t'+str(mutant_array[2])+'\t'+str(mutant_array[3]))
    else:
        print(str(mutant_array[0])+'\t'+str(mutant_array[1])+'\t'+str(mutant_array[2]))

# Load information outputed from previous cell and arranged in table manually
if estradiol == 100:
    back = pd.read_excel("./CoS_dgr_replicates_error_"+strain+"_full_estra.xlsx", index_col = 0)
elif estradiol == 0:
    back = pd.read_excel("./CoS_dgr_replicates_error_"+strain+"_full.xlsx", index_col = 0)
data_error = back.copy()
data_error[['MTX_dgr_max','MTX_dgr_min','DMSO_dgr_max','DMSO_dgr_min']] = np.nan
data_error = data_error.reset_index(drop = True)


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

# Construct plot with error bars
## Load data
source_plate = './'+strain+'_growthrate_vs_DMS.xlsx'
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


g2 = sns.scatterplot(data = df_wt, x = 'complementation', y='resistance', hue=['#d55e00',"#0072b2","#f0e442","#009e73"]
                     ,edgecolor="black", marker='X', s=100, legend=False, zorder=10,ax=ax)

norm = plt.Normalize(-0.75, 1.5)
sm = plt.cm.ScalarMappable(cmap="coolwarm",norm = norm)
cax = fig.add_axes([ax.get_position().x1+0.1, ax.get_position().y0+0.04, 0.02, ax.get_position().height])
ax.figure.colorbar(sm, cax=cax)

# Add other plot elements
coef, p = spearmanr(graph['complementation'], graph['resistance'])
coef = round(coef, 2)
if estradiol == 100:
    p = round(p, 15)
elif estradiol == 0:
    p = round(p, 14)
ax.text(s="ρ = {}\np = {}".format(coef, p), x=0.02, y=0.25, fontsize=20)
ax.set_ylabel("Growth rate - MTX IC90 (OD600$\cdot h^{-1}$)", fontsize=20)
ax.set_xlabel("Growth rate - DMSO (OD600$\cdot h^{-1}$)", fontsize=20)
ax.set_xlim(0,0.35)
ax.set_ylim(0,0.35)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', labelsize=20)
cax.tick_params(axis='y', labelsize=20)
ax.axline([ax.get_xlim()[0], ax.get_ylim()[0]], [ax.get_xlim()[1], ax.get_ylim()[1]], color = "black", linestyle = ":")
fig.tight_layout()
if estradiol == 100:
    plt.savefig("./FigureS8B.png", format = "png", dpi = 300, bbox_inches='tight')
elif estradiol == 0:
    plt.savefig("./FigureS8A.png", format = "png", dpi = 300, bbox_inches='tight')
    
    
#%%


# =============================================================================
# 
# #Similar concept but for data on IGA130 to compare growht rate between FDR and IGA validations
# 
# =============================================================================

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
# Define plot title order
order = pd.read_excel("./Validations_FDR0001_order.xlsx")
order = order['strain']

# Create copy of imported data
merged_df_tot = merged_df
# Set values within new dataframe to ease plotting
merged_df_tot['row'] = merged_df_tot['valid_id'].str[0]
merged_df_tot['col'] = merged_df_tot['valid_id'].str[1:]
merged_df_tot = merged_df_tot.sort_values(by=["row","col"])
#merged_df_tot = merged_df_tot.loc[~merged_df_tot['strain'].isin(["Empty","ScDHFR","PjDHFR",'mDHFR'])]
row_ord = ["A","B","C","D","E","F","G","H"]
col_ord = ["1","2","3","4","5","6"]

# Convert g/mL units into uM
# MTX molecular mass: 454.44g/mol
merged_df_tot['concentration_uM'] = merged_df_tot[merged_df_tot.drug=='MTX']['concentration'].apply(lambda x: convert_to_uM(x, 454.44))
merged_df_tot = merged_df_tot.fillna(0)
#%%
# Get growth rate from derivative

# Generate intermediate dataframe for data storage
temp_df = pd.DataFrame()
dg_df = pd.DataFrame(columns = ["well", "mutant", "dgr"])

# =============================================================================
concentration = 20 # 0 for DMSO, 20 for IC90
estra = 0 # 0 for no estra, 100 for estra
# =============================================================================

# Compute data to get growth rate and remove unecessary values
start_df = merged_df[merged_df["concentration"] == concentration]
#start_df = start_df[start_df["estradiol"] == estra]
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


## Growth curves were done with and without estradiol, specify conditions here ##
#=====================================================
estradiol = 0 # 0 for no estradiol, 100 for with estradiol
concentration = 20 # DMSO for no drug, 20 for IC90
#=====================================================
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

# Generate temporary variables
temp_df = pd.DataFrame()
dg_df = pd.DataFrame(columns = ["well", "mutant", "dgr"])
index = 0

## Set desired parameters
start_df = merged_df[merged_df["concentration"] == concentration]
#start_df = start_df[start_df["estradiol"] == estradiol]

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
    if len(mutant_array) > 3:
        print(str(mutant_array[0])+'\t'+str(mutant_array[1])+'\t'+str(mutant_array[2])+'\t'+str(mutant_array[3]))
    else:
        print(str(mutant_array[0])+'\t'+str(mutant_array[1])+'\t'+str(mutant_array[2]))

# Load information outputed from previous cell and arranged in table manually
error_FDR = pd.read_excel("./CoS_dgr_replicates_error_FDR0001_full.xlsx", index_col = 0)
error_IGA = pd.read_excel("./CoS_dgr_replicates_error_IGA130_full.xlsx", index_col = 0)
data_error_FDR = error_FDR.copy()
data_error_IGA = error_IGA.copy()
data_error_FDR[['MTX_dgr_max','MTX_dgr_min','DMSO_dgr_max','DMSO_dgr_min']] = np.nan
data_error_IGA[['MTX_dgr_max','MTX_dgr_min','DMSO_dgr_max','DMSO_dgr_min']] = np.nan
data_error_FDR = data_error_FDR.reset_index(drop = True)
data_error_IGA = data_error_IGA.reset_index(drop = True)


# Calculate data for error bars - FDR
medians = ['median_dgr-MTX', 'median_dgr-DMSO']
for row in data_error_FDR.index:
    data_mutant = data_error_FDR.iloc[row,]
    MTXdgr_median = data_mutant.at['median_dgr-MTX']
    array_dgr_MTX = [data_mutant.at['dgr1_MTX'],data_mutant.at['dgr2_MTX'],data_mutant.at['dgr3_MTX']]
    
    DMSOdgr_median = data_mutant.at['median_dgr-DMSO']
    array_dgr_DMSO = [data_mutant.at['dgr1_DMSO'],data_mutant.at['dgr2_DMSO'],data_mutant.at['dgr3_DMSO']]

    array_array = [array_dgr_MTX,array_dgr_DMSO]
    
    for array in array_array:
        if array == array_dgr_MTX:
            median = MTXdgr_median
            for value in array:
                if value > median:
                    up = value - median
                    data_error_FDR.at[row,'MTX_dgr_max'] = up
                if value < median:
                    down = median - value
                    data_error_FDR.at[row,'MTX_dgr_min'] = down
        if array == array_dgr_DMSO:
            median = DMSOdgr_median
            for value in array:
                if value > median:
                    up = value - median
                    data_error_FDR.at[row,'DMSO_dgr_max'] = up
                if value < median:
                    down = median - value
                    data_error_FDR.at[row,'DMSO_dgr_min'] = down
  
  
data_error_FDR = data_error_FDR.fillna(0)
data_error_FDR.drop(data_error_FDR.tail(1).index,inplace=True)

# Calculate data for error bars - IGA
medians = ['median_dgr-MTX', 'median_dgr-DMSO']
for row in data_error_IGA.index:
    data_mutant = data_error_IGA.iloc[row,]
    MTXdgr_median = data_mutant.at['median_dgr-MTX']
    array_dgr_MTX = [data_mutant.at['dgr1_MTX'],data_mutant.at['dgr2_MTX'],data_mutant.at['dgr3_MTX']]
    
    DMSOdgr_median = data_mutant.at['median_dgr-DMSO']
    array_dgr_DMSO = [data_mutant.at['dgr1_DMSO'],data_mutant.at['dgr2_DMSO'],data_mutant.at['dgr3_DMSO']]

    array_array = [array_dgr_MTX,array_dgr_DMSO]
    
    for array in array_array:
        if array == array_dgr_MTX:
            median = MTXdgr_median
            for value in array:
                if value > median:
                    up = value - median
                    data_error_IGA.at[row,'MTX_dgr_max'] = up
                if value < median:
                    down = median - value
                    data_error_IGA.at[row,'MTX_dgr_min'] = down
        if array == array_dgr_DMSO:
            median = DMSOdgr_median
            for value in array:
                if value > median:
                    up = value - median
                    data_error_IGA.at[row,'DMSO_dgr_max'] = up
                if value < median:
                    down = median - value
                    data_error_IGA.at[row,'DMSO_dgr_min'] = down
  
  
data_error_IGA = data_error_IGA.fillna(0)

# Construct plot with error bars
## Load data
source_plate = './FDR0001_vs_IGA130_growthrate.xlsx'
df = pd.read_excel(source_plate, header=0)
graph = pd.DataFrame()

# ------
if concentration == 0:
    condi = "DMSO"
elif concentration == 20:
    condi = "MTX"
# ------

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
graph["IGA130_MTX"] = resistance["IGA130"]
graph["FDR_MTX"] = resistance["FDR0001"]
graph["IGA130_DMSO"] = complementation["IGA130"]
graph["FDR_DMSO"] = complementation["FDR0001"]
graph["CoS_DMSO"] = complementation['CoS20']
graph["CoS_MTX"] = resistance['CoS20']
ctrls = ["Empty", "PjDHFR", 'ScDHFR']
df_wt = graph.loc[graph['mutant'].isin(ctrls)]
graph = graph.loc[~graph['mutant'].isin(ctrls)]

xa = "IGA130_" + condi
ya = "FDR_" + condi
cos = "CoS_" + condi

# Construct plot
fig, ax = plt.subplots(figsize = (8,8))
g1 = ax.scatter(data = graph, x = xa, y = ya, 
                    c=graph[cos], cmap="coolwarm", vmin=-0.75, vmax=1.5 
                    ,zorder=10,linewidth = 0)

minmin = condi + "_dgr_min"
xerror = [data_error_IGA[minmin].to_list(),data_error_IGA[minmin].to_list()]
yerror = [data_error_FDR[minmin].to_list(),data_error_FDR[minmin].to_list()]
ax.errorbar(x=data_error_IGA['median_dgr-'+condi], y=data_error_FDR['median_dgr-'+condi], xerr = xerror, yerr = yerror, color = 'black', fmt = 'o')

ctrl_hue = ['#d55e00',"#0072b2","#fff700"]
g2 = sns.scatterplot(data = df_wt, x = xa, y=ya, hue=ctrl_hue, palette=ctrl_hue, edgecolor="black", marker='X', s=100, legend=False, zorder=10,ax=ax)

norm = plt.Normalize(-0.75, 1.5)
sm = plt.cm.ScalarMappable(cmap="coolwarm",norm = norm)
cax = fig.add_axes([ax.get_position().x1+0.1, ax.get_position().y0+0.04, 0.02, ax.get_position().height])
ax.figure.colorbar(sm, cax=cax)

# Add other plot elements
coef, p = spearmanr(graph[xa], graph[ya])
coef = round(coef, 2)
if condi == "DMSO":
    p = round(p, 2)
elif condi == "MTX":
    p = round(p, 13)
ax.text(s="ρ = {}\np = {}".format(coef, p), x=0.02, y=0.25, fontsize=20)
ax.set_ylabel("Growth rate - FDR0001 (OD600$\cdot h^{-1}$)", fontsize=20)
ax.set_xlabel("Growth rate - IGA130 (OD600$\cdot h^{-1}$)", fontsize=20)
if condi == "DMSO":
    cax.set_ylabel("Selection coefficient - DMSO", fontsize=20, labelpad=15)
elif condi == "MTX":
    cax.set_ylabel("Selection coefficient - MTX IC90", fontsize=20, labelpad=15)
ax.set_xlim(0,0.35)
ax.set_ylim(0,0.35)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', labelsize=20)
cax.tick_params(axis='y', labelsize=20)
ax.axline([ax.get_xlim()[0], ax.get_ylim()[0]], [ax.get_xlim()[1], ax.get_ylim()[1]], color = "black", linestyle = ":")
#fig.tight_layout()
if condi == "DMSO":
    plt.savefig("./FigureS8C.png", format = "png", dpi = 300, bbox_inches='tight')
elif condi == "MTX":
    plt.savefig("./FigureS8D.png", format = "png", dpi = 300, bbox_inches='tight')
