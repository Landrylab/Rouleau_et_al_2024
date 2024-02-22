#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Francois D. Rouleau
Growth curves and dose response analysis, as well as pannels 1C and 1D.
Written on, and for, Spyder4.2.5

### Disclaimer
@Working directory is to be set in same location as script
@If working on something else than Spyder, #%% is what is used as cell separator. 
    Run each block of code/cell one after the other. 
"""
#Dose-response curves code
#%%
# Import relevant packages
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math
import scipy as sci
#%%
#Define important functions

##Get derivative growth rate
def get_derivative_growth_rate(g):
    # Rate is multiplied by 4 to convert from /15min (measurements every 15 min) to /h
    data_diff = g.OD.diff()*4
    # Get top 5 values 
    get_top_5 = data_diff.nlargest(5)    
    # Get median of top5
    growth_rate = np.median(get_top_5)
    return growth_rate

##Convert ug/mL to uM
def convert_to_uM(ug, mass_mol):
    # a simple function to convert ug/ml to uM with the molecular weigth of the compound
    vol = 0.001
    uM = (((ug*10**-6)/(vol*mass_mol))*10**6)
    return uM
#Metothrexate has a molecular mass of 454.44

##Get inhibition coefficient. Reference strain is Empty in no drug concentration
def get_inhib_coeff(df, fitness, fitness_type):
    # Inhibition coefficient is calculated either from the derivative growth rate (fitness_type = 'dgr')
    # or the normalized AUC parameter (fitness_type = 'norm_auc')
    ref = df.loc[(df.strain=='Empty') & (df.concentration_uM==0), fitness_type].mean()
    return (ref - fitness)/ref

##Calculate Hill equation
def hill_equation(x, E, n):
    # this function defines the hill equation, where:
    # x is the drug concentration
    # E is the concentration with half effect (IC50)
    # n is the hill coefficient,

    return 1/(1+((E/x)**n))

##Calculate IC concentrations
def get_IC(v, E, n):
    # This function calculates ICv based on IC50 (E) and hill coefficient (n) for any v value between 0 and 1
    print("Hill coefficient = "+str(round(n,2)))
    return E*((v/(1-v))**(1/n))


#%%
#Import plate plan
source_plate = './Plate_definition_growth_curves.xlsx'
plate_df = pd.read_excel(source_plate, header=0, index_col=0)
#%%
#Import Raw plate reader data 
source_table = './Raw_plate_read_data_growth_curves.xlsx'
source_df = pd.read_excel(source_table, header=1, index_col=0,skiprows=[2])
#Change seconds to hours
source_df.rename(columns=lambda x: float(x.strip('s'))/3600, inplace=True)

#%% 
#Put data in longform
longdf = source_df.reset_index().rename(columns = {'index' : 'well'})
longdf = longdf.melt(id_vars = ['well'], var_name = 'time (h)', value_name = 'OD')

#%%
#Annotate growthcurves using imported layout
merged_df = longdf.merge(right=plate_df, on=['well'])
#Change OD to log2
merged_df['log2_OD']= np.log2(merged_df['OD'])

#%%
#Generate figures with growth curves (Figure 1C)
#Change concentration to uM 
merged_df['concentration_uM'] = merged_df[merged_df.compound=='MTX']['concentration'].apply(lambda x: convert_to_uM(x, 454.44))
merged_df['concentration_uM'] = merged_df['concentration_uM'].round(2)

#Define plot data and structure
grid = sns.FacetGrid(data=merged_df[merged_df['time (h)'] <=40], col = 'cols',row='rows',col_order=[1,2],row_order=[1,2], hue='concentration_uM', palette='rocket_r',
                     height = 4.12 , aspect = 1.5)
grid.map(sns.lineplot, 'time (h)', 'OD')

#Define plot elements
grid.set_titles(col_template='{col_name}')
grid.set_axis_labels('Time (h)', '$\mathregular{OD_{600}}$', size = 20)
grid.tick_params(axis='both', which='major', labelsize=20)
grid.add_legend(title = ' MTX (μM)', fontsize = 20, loc='center right'#, borderaxespad = -1.2
                #, bbox_to_anchor=(0.6,0.5)
                , label_order = ["0.00","0.09","0.17","0.34","0.69","1.38","2.75","5.50","11.00","22.01","44.01","440.01"])
grid.fig.subplots_adjust(top=0.9)

#Set subplot titles
titles = ['Empty','PjDHFR','ScDHFR','mDHFR']
for ax, title in zip(grid.axes.flatten(), titles):
    ax.set_title(title, size = 20)
    
plt.setp(grid._legend.get_title(), fontsize=20)

grid.tight_layout()
#Save figure
plt.savefig('./Figure1C.png',format='png',dpi=200)
#%%

# =============================================================================
# 
# 
# Dose response curves and IC measurments
# 
# 
# =============================================================================

#%%
#Apply function to dataframe, create column with drg from top5 median
dgr = merged_df[merged_df['time (h)'] <= 40].groupby('well')[['OD']].apply(func=get_derivative_growth_rate).reset_index(name='dgr')
condensed_df = plate_df.merge(right=dgr, on=['well'])
#%%
#Set fitness type calculation and concentrations to uM
condensed_df['concentration_uM'] = condensed_df[condensed_df.compound=='MTX']['concentration'].apply(lambda x: convert_to_uM(x, 454.44))
condensed_df
fitness_type = 'dgr'
#%%
# Calculate inhibition coefficient and growth, correct for abberant values.
condensed_df['inhib_coeff'] = condensed_df[fitness_type].apply(lambda x: get_inhib_coeff(condensed_df, x, fitness_type))
#Get vales for growth (inverse of inhib)
condensed_df["growth"] = 1-condensed_df['inhib_coeff']
#Correct growth values greater than 1 to 1
condensed_df["growth"][condensed_df["growth"] > 1] = 1
#%%
#Generate fit data for the different strains used
fit_data_Sc_WT = condensed_df[(condensed_df.strain=='ScDHFR') & (condensed_df.concentration_uM!=0)].groupby(['concentration_uM'])[['inhib_coeff']].mean().reset_index()
fit_data_Murine = condensed_df[(condensed_df.strain=='mDHFR') & (condensed_df.concentration_uM!=0)].groupby(['concentration_uM'])[['inhib_coeff']].mean().reset_index()
fit_data_Pj_DHFR = condensed_df[(condensed_df.strain=='PjDHFR') & (condensed_df.concentration_uM!=0)].groupby(['concentration_uM'])[['inhib_coeff']].mean().reset_index()
fit_data_Empty = condensed_df[(condensed_df.strain=='Empty') & (condensed_df.concentration_uM!=0)].groupby(['concentration_uM'])[['inhib_coeff']].mean().reset_index()
#%%
# Calculate different values for curve fitting
wt_conc_WT = fit_data_Sc_WT.concentration_uM.values
wt_inhib_WT = fit_data_Sc_WT.inhib_coeff.values
wt_growth_WT = 1-wt_inhib_WT

wt_conc_Murine = fit_data_Murine.concentration_uM.values
wt_inhib_Murine = fit_data_Murine.inhib_coeff.values
wt_growth_Murine = 1-wt_inhib_Murine

wt_conc_Pj = fit_data_Pj_DHFR.concentration_uM.values
wt_inhib_Pj = fit_data_Pj_DHFR.inhib_coeff.values
wt_growth_Pj = 1-wt_inhib_Pj

wt_conc_Empty = fit_data_Empty.concentration_uM.values
wt_inhib_Empty = fit_data_Empty.inhib_coeff.values
wt_growth_Empty = 1-wt_inhib_Empty
#%%
#Fit Hill equation

#ScDHFR p0=[500,1]
# Get best fit parameters (popt) and covariance matrix (pcov)
popt_WT, pcov_WT = sci.optimize.curve_fit(hill_equation, wt_conc_WT, wt_growth_WT, p0=[500,1])
print(popt_WT)
print(pcov_WT)
#Murine p0=[2500,1]
popt_Murine, pcov_Murine = sci.optimize.curve_fit(hill_equation, wt_conc_Murine, wt_growth_Murine, p0=[2500,1])
print(popt_Murine)
print(pcov_Murine)
#Pj_DHFR p0=[2,1]
popt_Pj, pcov_Pj = sci.optimize.curve_fit(hill_equation, wt_conc_Pj, wt_growth_Pj, p0=[2,1])
print(popt_Pj)
print(pcov_Pj)

#Empty p0=[0.1,1]
popt_Empty, pcov_Empty = sci.optimize.curve_fit(hill_equation, wt_conc_Empty, wt_growth_Empty, p0=[0.1,1])
print(popt_Empty)
print(pcov_Empty)
#%%
# Calculate Rsquare and confidence intervals
# Use either: Empty, PjDHFR, ScDHFR, mDHFR
protein = "PjDHFR"
if protein == "PjDHFR":
    x = wt_conc_Pj
    y = wt_growth_Pj
    opt = popt_Pj
    cov = pcov_Pj
elif protein == "Empty":
    x = wt_conc_Empty
    y = wt_growth_Empty
    opt = popt_Empty
    cov = pcov_Empty
elif protein == "ScDHFR":
    x = wt_conc_WT
    y = wt_growth_WT
    opt = popt_WT
    cov = pcov_WT
elif protein == "mDHFR":
    x = wt_conc_Murine
    y = wt_growth_Murine
    opt = popt_Murine
    cov = pcov_Murine

# Calculate Rsquare
residuals = y- hill_equation(x, *opt)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((y-np.mean(y))**2)
r_squared = 1 - (ss_res / ss_tot)

print("Rsquared for "+protein+" is "+str(round(r_squared,3)))

# Compute IC confidence intervals
confidence = np.sqrt(np.diagonal(cov))
print("Uncertainty around ICs for "+protein+" is "+str(round(confidence[0]-confidence[1],2))+" uM")

#%%
#Get IC values for different strains and print them
ic10_WT = round(get_IC(0.9,popt_WT[0],popt_WT[1]),3)
print('IC10_WT = {} uM'.format(ic10_WT))
ic25_WT = round(get_IC(0.75,popt_WT[0],popt_WT[1]),3)
print('IC25_WT = {} uM'.format(ic25_WT))
ic50_WT = round(get_IC(0.5,popt_WT[0],popt_WT[1]),3)
print('IC50_WT = {} uM'.format(ic50_WT))
ic75_WT = round(get_IC(0.25,popt_WT[0],popt_WT[1]),3)
print('IC75_WT = {} uM'.format(ic75_WT))
ic10_Murine = round(get_IC(0.9,popt_Murine[0],popt_Murine[1]),3)
print('IC10_Murine = {} uM'.format(ic10_Murine))
ic25_Murine = round(get_IC(0.75,popt_Murine[0],popt_Murine[1]),3)
print('IC25_Murine = {} uM'.format(ic25_Murine))
ic50_Murine = round(get_IC(0.5,popt_Murine[0],popt_Murine[1]),3)
print('IC50_Murine = {} uM'.format(ic50_Murine))
ic75_Murine = round(get_IC(0.25,popt_Murine[0],popt_Murine[1]),3)
print('IC75_Murine = {} uM'.format(ic75_Murine))
ic10_Pj = round(get_IC(0.9,popt_Pj[0],popt_Pj[1]),3)
print('IC10_Pj = {} uM'.format(ic10_Pj))
ic25_Pj = round(get_IC(0.75,popt_Pj[0],popt_Pj[1]),3)
print('IC25_Pj = {} uM'.format(ic25_Pj))
ic50_Pj = round(get_IC(0.5,popt_Pj[0],popt_Pj[1]),3)
print('IC50_Pj = {} uM'.format(ic50_Pj))
ic75_Pj = round(get_IC(0.25,popt_Pj[0],popt_Pj[1]),3)
print('IC75_Pj = {} uM'.format(ic75_Pj))
ic90_Pj = round(get_IC(0.1,popt_Pj[0],popt_Pj[1]),3)
print('IC90_Pj = {} uM'.format(ic90_Pj))
ic10_Empty = round(get_IC(0.9,popt_Empty[0],popt_Empty[1]),3)
print('IC10_Empty = {} uM'.format(ic10_Empty))
ic25_Empty = round(get_IC(0.75,popt_Empty[0],popt_Empty[1]),3)
print('IC25_Empty = {} uM'.format(ic25_Empty))
ic50_Empty = round(get_IC(0.5,popt_Empty[0],popt_Empty[1]),3)
print('IC50_Empty = {} uM'.format(ic50_Empty))
ic75_Empty = round(get_IC(0.25,popt_Empty[0],popt_Empty[1]),3)
print('IC75_Empty = {} uM'.format(ic75_Empty))
ic90_Empty = round(get_IC(0.1,popt_Empty[0],popt_Empty[1]),3)
print('IC90_Empty = {} uM'.format(ic90_Empty))
#%%
#Generate values to make graph in log2
fit_vals_WT = np.logspace(np.log2(wt_conc_WT.min()),np.log2(wt_conc_WT.max()), num=50, base=2)
fit_vals_Murine = np.logspace(np.log2(wt_conc_Murine.min()),np.log2(wt_conc_Murine.max()), num=50, base=2)
fit_vals_Pj = np.logspace(np.log2(wt_conc_Pj.min()),np.log2(wt_conc_Pj.max()), num=50, base=2)
fit_vals_Empty = np.logspace(np.log2(wt_conc_Empty.min()),np.log2(wt_conc_Empty.max()), num=50, base=2)
#%%
#Generate Dose response curve figure

#Set general plot parameters
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

#Define plot area and content
fig, ax = plt.subplots(figsize=(14,8))
g1 = sns.scatterplot(data=condensed_df[condensed_df.concentration < 200], x='concentration_uM', y='growth', hue='strain', hue_order=["Empty", "PjDHFR", "ScDHFR", "mDHFR"], 
                     palette=['#d55e00','#0072b2','#f0e442','#009e73'], ax=ax)


#Define and plot fitted curves
ax.plot(fit_vals_Empty, [hill_equation(x, popt_Empty[0], popt_Empty[1]) for x in fit_vals_Empty], color='#d55e00')
ax.plot(fit_vals_Pj, [hill_equation(x, popt_Pj[0], popt_Pj[1]) for x in fit_vals_Pj], color='#0072b2')
ax.plot(fit_vals_WT, [hill_equation(x, popt_WT[0], popt_WT[1]) for x in fit_vals_WT], color='#f0e442')
ax.plot(fit_vals_Murine, [hill_equation(x, popt_Murine[0], popt_Murine[1]) for x in fit_vals_Murine], color='#009e73')

#Define plot elements
g1.set(ylim=(None, 1.05))
g1.set(xlim=(-1, 50))
g1.set_xlabel('Concentration (μM)', size = 22)
g1.set_ylabel('Growth coefficient', size = 22)
g1.legend(loc="center right", bbox_to_anchor=(1.035, 0.5), framealpha = 0, fontsize = 22, borderaxespad = 1,
          markerscale = 2, handletextpad = -0.2, title = "DHFR", title_fontsize = 22)  # Place legend outside the plot on center right
ax.text(10, 0.47, 'PjDHFR\nIC50 = '+str(round(ic50_Pj, 2))+' ± 0.23 μM\nN = -0.76', size = 22)
ax.text(25, 0.47, 'Empty\nIC50 = '+str(round(ic50_Empty, 2))+' ± 0.02 μM\nN = -0.50', size = 22)
ic90_Pj = round(get_IC(0.086,popt_Pj[0],popt_Pj[1]),2)
ax.tick_params(axis='both', which='major', labelsize=22)


# Save the figure
plt.tight_layout()
plt.savefig('./Figure1D.png', format='png', dpi=200)
