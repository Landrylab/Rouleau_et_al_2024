#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 10:51:52 2023

@author: francois
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.stats as stats
import seaborn as sns
import statsmodels
import sklearn.mixture
from upsetplot import from_indicators
from upsetplot import plot
#%%
### Ignore this, it is called in the next cell
longform = pd.read_csv("../Longform_dataframes/total_"+condition+"_longform_CoS_t15.csv")
#%%
# Welch t-test to find which samples are statistically significative
runcell(0)
aa_list = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", 
                "A", "V", "I", "L", "M", "F", "Y", "W", "*"]
seq ="DWQKSLTLIVALTLSRGIGLKNDLPWKLKSDMMFFSRVTSGLLVTRSTGQMNVVLMGRKTWESLPAHSRPLKNRINVVISRQEVLDLGGGAYHARSLDDALALLSQIYDSTSKIQLNRVFVIGGGELCKAAMEHSRLNRIIATVIHNEVDCDVFFPIDFRSSQSCLPWRKQDHSVLEAWVGSKVPQGKINENGFIYEFEMWIRD"
positions = list(range(2,206))

# Run cell to load condition data
conditions = ["MTX4","MTX20","DMSO"]
for condition in conditions:
    runcell(1)
    
    # Create silent dataframe as reference
    longform_silent = longform.copy()
    longform_silent = longform_silent[longform_silent['nature']=='silent']
    longform_silent = longform_silent.reset_index()
    cols_to_keep = [condition+'-1_CoS', condition+'-2_CoS', condition+'-3_CoS', 'index']
    cols_to_keep2 = [condition+'-1_CoS', condition+'-2_CoS', condition+'-3_CoS']
    
    # Rearrange dataframe and create new columns
    longform_silent_melt = longform_silent.copy()
    longform_silent_melt = longform_silent_melt[cols_to_keep]
    longform_silent_melt = longform_silent_melt.melt(id_vars=['index'], value_vars = cols_to_keep, ignore_index = False, var_name="Sample", value_name=condition+'_CoS')
    longform['p-values'] = 0.00
    longform['Welch_t-test'] = 0.00
    
    # Create silent array as reference and amino acid grouped array
    silent_array = longform_silent_melt[condition+'_CoS']
    longform_aa = pd.DataFrame(columns=['aa_WT','position','aa_mut','p-values','Welch_t-test','codons_num','median_CoS', 'mut_code', 'nature'])
    
    # Iterate through array to group amino acid mutants and do one sided greater Welch t-test
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

# Save data
# Only consider data with selection coefficient greater than 0 from Welch's t-test
    if condition == "MTX4":
        longform_aa_MTX4 = longform_aa.copy()
        greater_0 = longform_aa_MTX4[longform_aa_MTX4["median_CoS"] > 0]
    elif condition == "MTX20":
        longform_aa_MTX20 = longform_aa.copy()
        greater_0 = longform_aa_MTX20[longform_aa_MTX20["median_CoS"] > 0]


# Measure FDR using benjamini-hoeckber correction
    array = np.array(greater_0['p-values'].astype('float64'))
    reject, pvals_corrected, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(
        array, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

# Create column for values that passed the threshold
    greater_0['pass0.05'] = reject
    greater_0['pvalues_corr_0.05'] = pvals_corrected
# Create arrays for both conditions
    if condition == "MTX4":
        zeropointfive_ic75_bh = greater_0[greater_0['pass0.05'] == True]
    elif condition == "MTX20":
        zeropointfive_ic90_bh = greater_0[greater_0['pass0.05'] == True]

# Measure intersection of significance
intersect_bh = pd.merge(zeropointfive_ic75_bh, zeropointfive_ic90_bh, how ='inner', on =['mut_code'])

# Check which passed very resistant threshold as well
thresh = intersect_bh.copy()
thresh = thresh[thresh["median_CoS_x"] > 0.431753]
thresh = thresh[thresh["median_CoS_y"] > 0.549683]
#%%
# =============================================================================
# Gaussians
# =============================================================================
# Find number of dimensions which best recapitulate a give condition (Test AIC and BIC)

# Set condition ## This will generate all relevant information and figures at once
condition = "MTX4" #('MTX4' for IC75 and 'MTX20' for IC90)

# Set specific seed and test between 1 and 11 components
seed = 42

#Set data
if condition == "MTX4":
    longform = longform_aa_MTX4.copy()
elif condition == "MTX20":
    longform = longform_aa_MTX20.copy()

X = np.array(longform['median_CoS']).reshape(-1, 1)
N = np.arange(1, 11)
models = [None for i in range(len(N))]
for i in range(len(N)):
    models[i] = sklearn.mixture.GaussianMixture(n_components=N[i], random_state = seed).fit(X, y=None)
# compute the AIC and the BIC
AIC = [m.aic(X) for m in models]
BIC = [m.bic(X) for m in models]

# Define plot information
font = {'family' : 'normal',
        'size'   : 20}
matplotlib.rc('font', **font)

# Create plot
fig = plt.figure(figsize=(6, 6))
fig.subplots_adjust(left=0.12, right=0.97,
                    bottom=0.21, top=0.9, wspace=0.5)
ax = fig.add_subplot(111)
ax.plot(N, AIC, '-k', label='AIC')
ax.plot(N, BIC, '--k', label='BIC')
ax.set_xlabel('Number of components')
ax.set_ylabel('Information criterion')
ax.legend(loc=1)

if condition == "MTX20":
    fig.savefig('./Fig_S1B.png', format = 'png', dpi = 200, bbox_inches='tight')
elif condition == "MTX4":
    fig.savefig('./Fig_S1A.png', format = 'png', dpi = 200, bbox_inches='tight')
plt.show()


# Once you have the number of components, get relevant information from thresholds
# Run sklearn to estimate Gaussian mixture
print(condition)
if condition == "MTX20":
    X = np.array(longform_aa_MTX20['median_CoS']).reshape(-1, 1)
    # Set specific seed for reproductibilty
    seed = 42
    
    # Run sklearn
    models3 = sklearn.mixture.GaussianMixture(n_components=4, random_state = seed).fit(X, y=None)
    test3 = models3.predict(X)
    
    # Return gaussian info to initial dataframe
    longform_aa_MTX20['gaussian'] = test3
    longform_aa_gaussian_mid = longform_aa_MTX20[longform_aa_MTX20['gaussian'] == 3]
    longform_aa_gaussian_top = longform_aa_MTX20[longform_aa_MTX20['gaussian'] == 1]
    
    # Get different info from micture model
    pose = list(dict.fromkeys(longform_aa_gaussian_top['position']))
    median_mid = longform_aa_gaussian_mid['median_CoS'].median()
    print(median_mid)
    pass_mid = longform_aa[longform_aa['median_CoS']>median_mid]
    print(len(pass_mid))
    median_top = longform_aa_gaussian_top['median_CoS'].min()
    print(median_top)
    pass_top = longform_aa[longform_aa['median_CoS']>median_top]
    print(len(pass_top))
    pose_mid_20 = list(dict.fromkeys(longform_aa_gaussian_mid['position']))
    pose_top_20 = list(dict.fromkeys(longform_aa_gaussian_top['position']))

elif condition == "MTX4":
    X = np.array(longform_aa_MTX4['median_CoS']).reshape(-1, 1)
    # Set specific seed for reproductibilty
    seed = 42
    
    # Run sklearn
    models3 = sklearn.mixture.GaussianMixture(n_components=5, random_state = seed).fit(X, y=None)
    test3 = models3.predict(X)
    
    # Return gaussian info to initial dataframe
    longform_aa['gaussian'] = test3
    longform_aa_gaussian_mid = longform_aa[longform_aa['gaussian'] == 1]
    longform_aa_gaussian_top = longform_aa[longform_aa['gaussian'] == 3]
    
    # Get different info from micture model
    median = longform_aa_gaussian_mid['median_CoS'].median()
    print(median)
    pass_mid = longform_aa[longform_aa['median_CoS']>median]
    print(len(pass_mid))
    median_top = longform_aa_gaussian_top['median_CoS'].min()
    print(median_top)
    pass_top = longform_aa[longform_aa['median_CoS']>median_top]
    print(len(pass_top))
    pose_mid_4 = list(dict.fromkeys(longform_aa_gaussian_mid['position']))
    pose_top_4 = list(dict.fromkeys(longform_aa_gaussian_top['position']))


# =============================================================================
# Generate plots for Gaussian mixture models
# =============================================================================
# Gaussian with Density mapping for BH  corrected

# Define condition
#Define plotting properties
font = {'family' : 'normal',
        'size'   : 20}
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)
matplotlib.rc('font', **font)

# Define othe parameters for plotting
if condition == "MTX20":
    mods = 3
    mini = -0.570587
    maxi = 1.63917
    passed_bh = zeropointfive_ic90_bh['median_CoS']
    ic = "IC90"
    pan = "S3D"
    longform = pd.read_csv("../Longform_dataframes/total_MTX20_longform_CoS_t15.csv")
elif condition == "MTX4":
    mods = 4
    mini = -0.943674
    maxi = 1.13987
    passed_bh = zeropointfive_ic75_bh['median_CoS']
    ic = "IC75"
    pan = "S3C"
    longform = pd.read_csv("../Longform_dataframes/total_MTX4_longform_CoS_t15.csv")

fig, ax = plt.subplots(figsize = (6,6), layout = 'tight')
M_best = models[mods]
#plt.suptitle('Best mixture model \n'+ic+' - '+str(mods+1)+' components', size =15)
x = np.linspace(mini, maxi, len(longform_aa_MTX20))
logprob = M_best.score_samples(x.reshape(-1, 1))
responsibilities = M_best.predict_proba(x.reshape(-1, 1))
pdf = np.exp(logprob)
pdf_individual = responsibilities * pdf[:, np.newaxis]

#Map gaussians
ax.hist(X, 30, density=True, histtype='stepfilled', alpha=0.4)
ax.plot(x, pdf, '-k')
ax.plot(x, pdf_individual, '--k')
ax.set_xlabel('$x$')
ax.set_ylabel('$p(x)$')
ax.set_xlim(-1,1.7)

#Map BH
ax1 = ax.twinx()
bh = sns.kdeplot(passed_bh, color = "green", ax=ax1, cut = 0)
ax1.set_ylim(0,len(longform)/len(passed_bh))
ax1.set(yticklabels=[])
ax1.set(ylabel=None)
ax1.tick_params(right=False)

ax.set_xlabel("Selection coefficient - MTX "+ic)
ax.set_ylabel("Density")
# =============================================================================
# if condition == "MTX20":
#     fig.savefig('./Fig_S1D.png', format = 'png', dpi = 200, bbox_inches='tight')
# elif condition == "MTX4":
#     fig.savefig('./Fig_S1C.png', format = 'png', dpi = 200, bbox_inches='tight')
# =============================================================================
plt.show()

#%%
# Upset plot to check intersection between the two conditions after FDR

# Create new dataframes with relevant info and remove duplicates
upset_df = pd.DataFrame()
ic75_bh_mut = zeropointfive_ic75_bh['mut_code']
ic90_bh_mut = zeropointfive_ic90_bh['mut_code']
all_muts = pd.concat([ic75_bh_mut,ic90_bh_mut])
all_muts = all_muts.drop_duplicates()
upset_df["Mutants"] = all_muts

# Format data to make it compatible with upset plot package
ic75_bh_bool = all_muts.isin(ic75_bh_mut)
ic90_bh_bool = all_muts.isin(ic90_bh_mut)
upset_df["IC75_BH"] = ic75_bh_bool
upset_df["IC90_BH"] = ic90_bh_bool
upset_df_cp = upset_df.copy()
upset_df_cp.set_index('Mutants', drop = True, inplace = True)
upset_df_form = from_indicators(upset_df_cp)

# Generate simple upset plot
plot(upset_df_form,subset_size="count")
plt.tight_layout()
plt.savefig('./Fig_S6.png', format = "png", dpi = 300)