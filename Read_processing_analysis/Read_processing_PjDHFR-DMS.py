#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 09:53:19 2023

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
import re
from Bio import SeqIO
from subprocess import call
import math
import matplotlib.patches as mpatches
import statistics
#%%
# Data was processed in a per-fragment manner. Repeat this analysis per fragment as necessary
os.chdir("Your directory")

# Import file list. Make manually using all of your relevant files
fragment = "F3" #F1, F2 or F3
os.chdir("./"+fragment)
R1=pd.read_csv('./'+fragment+'_list_R1.csv', header = None)
file_list1 = R1.loc[0, :].values.tolist()
R2=pd.read_csv('./'+fragment+'_list_R2.csv', header = None)
file_list2 = R2.loc[0, :].values.tolist()
#%%
# Run fastqc on all files, make sure read quality is ok
for file in file_list1:
    path = "./"+file+".fastq"
    transeq_call = "fastqc "+path+ " --outdir=./fastqc_"+fragment+"/"
    subprocess.check_output(transeq_call, shell=True)
     
for file in file_list2:
    path = "./"+file+".fastq"
    transeq_call = "fastqc "+path+ " --outdir=./fastqc_"+fragment+"/"
    subprocess.check_output(transeq_call, shell=True)


#%%
# Merge reads using pandaseq ## This version of the command also uses
# cutadapt-like parameters to trim adapters using alignments ## 

## Create directory for trimmed/merged reads
if not os.path.exists("./merge_pre-trim_"+fragment):
    os.makedirs("./merge_pre-trim_"+fragment)

## Set adapter sequence to be removed to only keep relevant data 
if fragment == "F1" :
    forw = "aacaccagaacttagtttcgacgg"
    rev = "TGGAGATGACAACGTTGATTCTGTT"
elif fragment == "F2" :
    forw = "CTTTACCAGCTCATTCTAGACCATTAAAG"
    rev = "GACATCACAGTCAACTTCGTTGTGAAT"
elif fragment == "F3" :
    forw = "GATTGAACCGTATCATTGCCACCGTC"
    rev = "tgactcgaggtcgacggtatcgattta"
else:
    pass

# Run pandaseq call on all R1 and R2 files of interest 
index = 0
for file in file_list1:
    pathR1 = "./"+file+".fastq"
    pathR2 = "./"+file_list2[index]+".fastq"
    out = "./merge_pre-trim_"+fragment+"/"+file+"_merged.fastq"
    panda_seq_call = 'pandaseq -f '+pathR1+' -r '+pathR2+ ' -L 550 -o 75 -O 400 -k 2 -B -N -t 0.6 -T 8 -p '+forw+' -q '+rev+' -w '+ out
    subprocess.check_output(panda_seq_call, shell=True)
    index = index + 1

#%%
# Aggregate trimmed and merged reads using vsearch for count
if not os.path.exists("./merge_pre-trim_"+fragment+"/aggregate"):
    os.makedirs("./merge_pre-trim_"+fragment+"/aggregate")

# Run vsearch call on all files of interest 
for file in file_list1:
    path = "./merge_pre-trim_"+fragment+"/"+file+"_merged.fastq"
    aggregate_call = "vsearch --derep_fulllength "+path+" --relabel seq --output "+"./merge_pre-trim_"+fragment+"/aggregate/"+file+"_agg.fasta"+" --sizeout"
    subprocess.check_output(aggregate_call, shell=True)

#%%
# Get file metadata
metadata = pd.DataFrame(columns = ["sample", "total", "unique", "singleton", "non-singleton", "fragment","replicate","condition"])
index = 0
sample_list = list()
unique_list = list()

for i in file_list1:
    #Find samplename
    sample = i.split("_")[0]
    sample_list.append(sample)
    metadata.loc[index,"sample"] = str(sample)
    #Find replicate id
    rep = i.split("-")[1]
    metadata.loc[index,"replicate"] = int(rep)
    #Find condition 
    cond1 = i.split("_S")[0]
    cond = cond1.split("-")[2]
    metadata.loc[index,"condition"] = str(cond)
    #Number of unique reads
    cmd_unique = "grep -c seq ./merge_pre-trim_"+fragment+"/aggregate/"+i+ "_agg.fasta"
    unique_num = subprocess.check_output(cmd_unique, shell=True)
    unique_num = int(re.search(r'\d+', str(unique_num)).group())
    unique_list.append(int(unique_num))
    metadata.loc[index,'unique'] = int(unique_num)
    metadata.loc[index, "fragment"] = i.split("-")[0]
    index = index + 1
    
    
    
    
# Get metadata for singleton reads/mutants to evaluate fraction of non-usable reads
total = 0
single = 0
index = 0

for i in file_list1:
    with open("./merge_pre-trim_"+fragment+"/aggregate/"+i+"_agg.fasta", 'r') as source:
        for line in source:
            if line.startswith('>')==True:
                seq_info = line.split(';')
                seq_id = seq_info[0]
                seq_count = int(line.split("=")[1])
                if int(seq_count) == 1:
                    single = single + 1
                    total = total + 1
                else:
                    total = total + seq_count
                    #print(list_single)
    metadata.loc[index,'singleton'] = int(single)
    metadata.loc[index, "total"] = int(total)
    metadata.loc[index, "non-singleton"] = metadata.loc[index, "unique"] - int(single)
    total = 0
    single = 0
    index = index+1
    

metadata.to_csv("../"+fragment+"_metadata_pandaseqTrim.csv")


    
#%%
# Translate Nt to aa using transeq
if not os.path.exists("./merge_pre-trim_"+fragment+"/aggregate/aa"):
    os.makedirs("./merge_pre-trim_"+fragment+"/aggregate/aa")
    

for file in file_list1:
    path = "./merge_pre-trim_"+fragment+"/aggregate/"+file+"_agg.fasta"
    transeq_call = "transeq "+path+" ./merge_pre-trim_"+fragment+"/aggregate/aa/"+file+"_aa.fasta -frame=1 -stdout Y"
    subprocess.check_output(transeq_call, shell=True)

#%%

# =============================================================================
# Define functions to be used later
# =============================================================================


def fasta2dict(fil):
    
#    Read fasta-format file, return dict of form scaffold:sequence.
#    Note: Uses only the unique identifier of each sequence, rather than the 
#    entire header, for dict keys. 
    
    dic = {}
    cur_scaf = ''
    cur_seq = []
    for line in open(fil):
        if line.startswith(">") and cur_scaf == '':
            cur_scaf = line.split('\n')[0]
        elif line.startswith(">") and cur_scaf != '':
            dic[cur_scaf] = ''.join(cur_seq)
            cur_scaf = line.split('\n')[0]
            cur_seq = []
        else:
            cur_seq.append(line.rstrip())
    dic[cur_scaf] = ''.join(cur_seq)
    return dic

def get_key(val, dicti):
    for key, value in dicti.items():
        if val == value:
            return key
        else:
            pass


def translate(seq):
      
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    protein =""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein+= table[codon]
    return protein

#%% 

## Process reads and create output file

### Step 1: Create dataframe with all sequences found in T0 and their frequency
df_T0 = pd.DataFrame()

#### 1.1: Create a list of all possible sequences found in T0
list_T0_tot = list()
T0_files = [s for s in file_list1 if "T0" in s]
for file in T0_files:
    dicti = fasta2dict("./merge_pre-trim_"+fragment+"/aggregate/"+file+"_agg.fasta")
    my_list = list(dicti.values())
    list_T0_tot.extend(my_list)

#### 1.2: Use dictionnary to dereplicate this list and add it to dataframe
list_T0_tot = list(dict.fromkeys(list_T0_tot))
df_T0["nt_seq"] = list_T0_tot
df_T0["aa_seq"] = ""

#### 1.3: Count frequency in all samples
for file in T0_files:
    rep = file.split("-")[1]
    samp = "T0-"+rep
    df_T0[samp] = ""
    dicti = fasta2dict("./merge_pre-trim_"+fragment+"/aggregate/"+file+"_agg.fasta")
    for line in list(range(len(df_T0))):
        seq = df_T0.at[line,"nt_seq"]
        key = get_key(seq, dicti)
        if key is not None:
            count = key.split("=")[1]
            df_T0.at[line,samp] = count
        print(line)
   
    
for line in list(range(len(df_T0))):
    seq = df_T0.at[line,"nt_seq"]
    prot = translate(seq)
    df_T0.at[line,"aa_seq"] = prot
    

#%%

# Create individual dataframes for all conditions
# In this code, MTX conditions are presented in ug/mL. MTX4 (ug/mL) is IC75, and MTX20 (ug/mL) is IC90.

df_T0_D = df_T0.copy()
df_T0_M4 = df_T0.copy()
df_T0_M20 = df_T0.copy()

conditions = ["DMSO", "MTX4", "MTX20"]
directory = "./merge_pre-trim_"+fragment+"/aggregate/"

# For all conditions, create a dataframe with all the possible T0 sequences

for cond in conditions:
    matching = [s for s in file_list1 if cond in s]
    if cond == "DMSO":
        df = df_T0_D
    elif cond == "MTX4":
        df = df_T0_M4
    elif cond == "MTX20":
        df = df_T0_M20
    for fil in matching:
        rep = fil.split("-")[1]
        cond_rep = cond+"-"+rep
        df[cond_rep] = ""
        dicti = fasta2dict("./merge_pre-trim_"+fragment+"/aggregate/"+fil+"_agg.fasta")
        for line in list(range(len(df))):
            seq = df.at[line,"nt_seq"]
            key = get_key(seq, dicti)
            print(line)
            if key is not None:
                count = key.split("=")[1]
                df.at[line,cond_rep] = int(count)
            else:
                count = 0
                df.at[line,cond_rep] = count

    if cond == "DMSO":
        df_T0_D = df
    elif cond == "MTX4":
        df_T0_M4 = df
    elif cond == "MTX20":
        df_T0_M20 = df

    print(cond + " is done")


df_T0_D_back = df_T0_D.copy()
df_T0_M4_back = df_T0_M4.copy()
df_T0_M20_back = df_T0_M20.copy()


#%%
### Do the mutation calling from dictionary

### Before doing this, run Generate_all_possible_mutants.py. This will generate all possible
### mutant for a given fragment. Following code will match assembled reads to all possible.
#### This part of the code is what does the sequence matching to mutant
all_possible = fasta2dict("../../"+fragment+"_all_variants.fasta")

for cond in conditions:
    matching = [s for s in file_list1 if cond in s]
    if cond == "DMSO":
        df = df_T0_D
    elif cond == "MTX4":
        df = df_T0_M4
    elif cond == "MTX20":
        df = df_T0_M20
    
    df["codon_WT"] = ""
    df["aa_WT"] = ""
    df["position"] = ""
    df["aa_mut"] = ""
    df["codon_mut"] = ""

    lines = list(range(len(df.index)))
    for line in lines:
        print(line)
        seq = df.iloc[line,0]
        mut_tot = get_key(seq, all_possible)
        if mut_tot is not None:
            mut_info = mut_tot[1:]
            mut_info_list = mut_info.split(";")
    
            df.at[line, "codon_WT"] = str(mut_info_list[0])
            df.at[line, "aa_WT"] = str(mut_info_list[1])
            df.at[line, "position"] = int(mut_info_list[2])
            df.at[line, "aa_mut"] = str(mut_info_list[3])
            df.at[line, "codon_mut"] = str(mut_info_list[4])


    if cond == "DMSO":
        df_T0_D = df
    elif cond == "MTX4":
        df_T0_M4 = df
    elif cond == "MTX20":
        df_T0_M20 = df

    print(cond + " is done")

# Make backups
df_T0_D_back = df_T0_D.copy()
df_T0_M4_back = df_T0_M4.copy()
df_T0_M20_back = df_T0_M20.copy()
#%%
# Reset DF to backup
df_T0_D = df_T0_D_back.copy()
df_T0_M4 = df_T0_M4_back.copy()
df_T0_M20 = df_T0_M20_back.copy()

#%%
# Savefiles as all_detected
df_T0_D.to_csv("./"+fragment+"_DMSO_all_detected.csv")
df_T0_M4.to_csv("./"+fragment+"_MTX4_all_detected.csv")
df_T0_M20.to_csv("./"+fragment+"_MTX20_all_detected.csv")


#%%
# Go through all mutant to make sure that not sequences with more than one mutated codon
# make it through to the alignement part 

fragment = "FX"
conditions = ["DMSO", "MTX4", "MTX20"]
df_T0_D = pd.read_csv("./"+fragment+"_DMSO_all_detected.csv", index_col=0)
df_T0_M4 = pd.read_csv("./"+fragment+"_MTX4_all_detected.csv", index_col=0)
df_T0_M20= pd.read_csv("./"+fragment+"_MTX20_all_detected.csv", index_col=0)
# Drop all values with aa hamming >1
for cond in conditions:
    matching = [s for s in file_list1 if cond in s]
    if cond == "DMSO":
        df = df_T0_D
    elif cond == "MTX4":
        df = df_T0_M4
    elif cond == "MTX20":
        df = df_T0_M20
    df = df.replace(r'^\s*$', np.nan, regex=True)


    len_pre = len(df.index)
    df.dropna(subset = ["codon_WT"], inplace = True)
    len_post = len(df.index)
    print((str(len_pre-len_post)+" lines have been dropped because of NA (hamming > 3 nt)"))
    df = df.reset_index(drop = True)
    if cond == "DMSO":
        df_T0_D = df
    elif cond == "MTX4":
        df_T0_M4 = df
    elif cond == "MTX20":
        df_T0_M20 = df


#%%
# Measure log2 fold changes:
repli = [1,2,3]
conditions = ["DMSO", "MTX4", "MTX20"]

# Replace all 0 with 1 to allow for Log2 calculation
# Calculate frequencies
# ====================================
# For T0
for cond in conditions:
    if cond == "DMSO":
        df = df_T0_D
    elif cond == "MTX4":
        df = df_T0_M4
    elif cond == "MTX20":
        df = df_T0_M20
    df.fillna(0, inplace = True)
    df=df.replace(0,1)
    for rep in repli:
        sam = "T0-"+str(rep)
        sam_freq = sam+"_freq"
        tot = df[sam].sum()
        for li in list(range(len(df))):
            rc = df.at[li, sam]
            freq = rc / tot 
            df.at[li, sam_freq] = freq
    if cond == "DMSO":
        df_T0_D = df
    elif cond == "MTX4":
        df_T0_M4 = df
    elif cond == "MTX20":
        df_T0_M20 = df
        

# For selection
for cond in conditions:
    if cond == "DMSO":
        df = df_T0_D
    elif cond == "MTX4":
        df = df_T0_M4
    elif cond == "MTX20":
        df = df_T0_M20
    

    sel1_total = df[cond+"-1"].sum()
    sel2_total = df[cond+"-2"].sum()
    sel3_total = df[cond+"-3"].sum()
    
    for rep in repli:
        sam = str(cond)+"-"+str(rep)
        sam_freq = sam+"_freq"
        tot = df[sam].sum()
        df[sam_freq] = ""
        
        for li in list(range(len(df))):
            rc = df.at[li, sam]
            freq = rc / tot 
            df.at[li, sam_freq] = float(freq)

    if cond == "DMSO":
        df_T0_D = df
    elif cond == "MTX4":
        df_T0_M4 = df
    elif cond == "MTX20":
        df_T0_M20 = df



# ====================================
# Calculate Log2freq
for cond in conditions:
    if cond == "DMSO":
        df = df_T0_D
    elif cond == "MTX4":
        df = df_T0_M4
    elif cond == "MTX20":
        df = df_T0_M20
    matching = ["T0-1_freq", "T0-2_freq", "T0-3_freq", cond+"-1_freq", cond+"-2_freq", cond+"-3_freq"]

    for ro in matching:
        sam_log = ro+"_log2"
        df[sam_log] = ""
        
        for li in list(range(len(df))):
            frequ = df.at[li, ro]
            log = math.log2(frequ)
            df.at[li, sam_log] = float(log)
    
    if cond == "DMSO":
        df = df_T0_D
    elif cond == "MTX4":
        df = df_T0_M4
    elif cond == "MTX20":
        df = df_T0_M20

# ====================================
#Calculate log2FC per sample
for cond in conditions:
    if cond == "DMSO":
        df = df_T0_D
    elif cond == "MTX4":
        df = df_T0_M4
    elif cond == "MTX20":
        df = df_T0_M20
        
    for rep in repli:
        pre = "T0-"+str(rep)+"_freq_log2"
        pos = cond+"-"+str(rep)+"_freq_log2"
        fc = cond+"-"+str(rep)+"_log2FC"
        df[fc] = ""
        
        for li in list(range(len(df))):
            foldc = df.at[li, pos] - df.at[li, pre]
            df.at[li, fc] = float(foldc)
    
    if cond == "DMSO":
        df = df_T0_D
    elif cond == "MTX4":
        df = df_T0_M4
    elif cond == "MTX20":
        df = df_T0_M20

# ====================================
# Calculate median 
for cond in conditions:
    if cond == "DMSO":
        df = df_T0_D
    elif cond == "MTX4":
        df = df_T0_M4
    elif cond == "MTX20":
        df = df_T0_M20
        
    rep1 = cond+"-1_log2FC"
    rep2 = cond+"-2_log2FC"
    rep3 = cond+"-3_log2FC"
    meda = "median_FC"
    df[meda] = ""
    
    for li in list(range(len(df))):
        media = pd.Series([df.at[li,rep1],df.at[li,rep2],df.at[li,rep3]])
        median = media.median()
        df.at[li, meda] = float(median)
    
    if cond == "DMSO":
        df = df_T0_D
    elif cond == "MTX4":
        df = df_T0_M4
    elif cond == "MTX20":
        df = df_T0_M20

#%%
# Save data to new file
df_T0_D.to_csv("./"+fragment+"_DMSO_all_matched.csv")
df_T0_M4.to_csv("./"+fragment+"_MTX4_all_matched.csv")
df_T0_M20.to_csv("./"+fragment+"_MTX20_all_matched.csv")


#%%
### Calculate selection coefficient

fragment = "total"
condition = "MTX20"


#Add number of generations by samples
longform_F1 = pd.read_csv("./F1/F1_"+condition+"_all_matched.csv", index_col=0)
longform_F1["fragment"] = "F1"
if condition == "DMSO":
    longform_F1['DMSO-1_gen'] = float(5.730096466)
    longform_F1['DMSO-2_gen'] = float(5.526068812)
    longform_F1['DMSO-3_gen'] = float(5.615004242)
elif condition == "MTX4":
    longform_F1['MTX4-1_gen'] = float(5.455162801)
    longform_F1['MTX4-2_gen'] = float(5.510961919)
    longform_F1['MTX4-3_gen'] = float(5.444600814)
elif condition == "MTX20":
    longform_F1['MTX20-1_gen'] = float(5.468909092)
    longform_F1['MTX20-2_gen'] = float(5.496334513)
    longform_F1['MTX20-3_gen'] = float(5.408032466)
else:
    pass
longform_F2 = pd.read_csv("./F2/F2_"+condition+"_all_matched.csv", index_col=0)
longform_F2["fragment"] = "F2"
if condition == "DMSO":
    longform_F2['DMSO-1_gen'] = float(5.696272084)
    longform_F2['DMSO-2_gen'] = float(5.659639187)
    longform_F2['DMSO-3_gen'] = float(5.6043679)
elif condition == "MTX4":
    longform_F2['MTX4-1_gen'] = float(5.350850805)
    longform_F2['MTX4-2_gen'] = float(5.247927513)
    longform_F2['MTX4-3_gen'] = float(5.240314329)
elif condition == "MTX20":
    longform_F2['MTX20-1_gen'] = float(5.371907298)
    longform_F2['MTX20-2_gen'] = float(5.394033895)
    longform_F2['MTX20-3_gen'] = float(5.269781238)
else:
    pass
longform_F3 = pd.read_csv("./F3/F3_"+condition+"_all_matched.csv", index_col=0)
longform_F3["fragment"] = "F3"
if condition == "DMSO":
    longform_F3['DMSO-1_gen'] = float(5.377123749)
    longform_F3['DMSO-2_gen'] = float(5.231509211)
    longform_F3['DMSO-3_gen'] = float(5.262658655)
elif condition == "MTX4":
    longform_F3['MTX4-1_gen'] = float(5.325890061)
    longform_F3['MTX4-2_gen'] = float(5.330199833)
    longform_F3['MTX4-3_gen'] = float(5.241458871)
elif condition == "MTX20":
    longform_F3['MTX20-1_gen'] = float(5.524189078)
    longform_F3['MTX20-2_gen'] = float(5.546585829)
    longform_F3['MTX20-3_gen'] = float(5.533874777)
else:
    pass


longform_F1["median_CoS"] = ""
longform_F1[condition+"-1_CoS"] = ""
longform_F1[condition+"-2_CoS"] = ""
longform_F1[condition+"-3_CoS"] = ""
longform_F1["fragment"] = "F1"
longform_F2["median_CoS"] = ""
longform_F2[condition+"-1_CoS"] = ""
longform_F2[condition+"-2_CoS"] = ""
longform_F2[condition+"-3_CoS"] = ""
longform_F2["fragment"] = "F2"
longform_F3["median_CoS"] = ""
longform_F3[condition+"-1_CoS"] = ""
longform_F3[condition+"-2_CoS"] = ""
longform_F3[condition+"-3_CoS"] = ""
longform_F3["fragment"] = "F3"


longform_int = longform_F1.append(longform_F2)
longform = longform_int.append(longform_F3)
#longform = longform[longform["position"] > 1]
longform = longform[longform["position"] < 206]
longform = longform.reset_index(drop = True)



aa_list = ["G","A","L","M","F","W","K","Q","E","S",
           "P","V","I","C","Y","H","R","N","D","T","*"]

aa_NNK_matrix = [
              'A-GCG', 'A-GCT',
              'C-TGT',
              'D-GAT',
              'E-GAG',
              'F-TTT',
              'G-GGG', 'G-GGT',
              'H-CAT',
              'I-ATT',
              'K-AAG',
              'L-CTG', 'L-CTT', 'L-TTG',
              'M-ATG',
              'N-AAT',
              'P-CCG', 'P-CCT',
              'Q-CAG',
              'R-AGG', 'R-CGG', 'R-CGT',
              'S-AGT', 'S-TCG', 'S-TCT',
              'T-ACG', 'T-ACT', 
              'V-GTG', 'V-GTT',
              'W-TGG',
              'Y-TAT',
              '*-TAG']
NNK_aa_dict = {
            'GCG':'A','GCT':'A',
            'TGT':'C',
            'GAT':'D',
            'GAG':'E',
            'TTT':'F',
            'GGG':'G','GGT':'G',
            'CAT':'H',
            'ATT':'I',
            'AAG':'K',
            'CTG':'L','CTT':'L','TTG':'L',
            'ATG':'M',
            'AAT':'N',
            'CCG':'P','CCT':'P',
            'CAG':'Q',
            'AGG':'R','CGG':'R','CGT':'R',
            'AGT':'S','TCG':'S','TCT':'S',
            'ACG':'T','ACT':'T',
            'GTG':'V','GTT':'V',
            'TGG':'W',
            'TAT':'Y',
            'TAG':'*'}


NNK_matrix = [
              'GCG', 'GCT',
              'TGT',
              'GAT',
              'GAG',
              'TTT',
              'GGG', 'GGT',
              'CAT',
              'ATT',
              'AAG',
              'CTG', 'CTT', 'TTG',
              'ATG',
              'AAT',
              'CCG', 'CCT',
              'CAG',
              'AGG', 'CGG', 'CGT',
              'AGT', 'TCG', 'TCT',
              'ACG', 'ACT',
              'GTG', 'GTT',
              'TGG',
              'TAT',
              'TAG']


if fragment == "F1" :
    lower = 2
    upper = 74
elif fragment == "F2" :
    lower = 74
    upper = 146
elif fragment == "F3" :
    lower = 146
    upper = 206
elif fragment == "total" :
    lower = 2
    upper = 206
else:
    pass


#Compound codons into aa (calculate median for codons)
if fragment == "F1" :
    seq = "DWQKSLTLIVALTLSRGIGLKNDLPWKLKSDMMFFSRVTSGLLVTRSTGQMNVVLMGRKTWESLPAHSRPLK"
elif fragment == "F2" :
    seq = "NRINVVISRQEVLDLGGGAYHARSLDDALALLSQIYDSTSKIQLNRVFVIGGGELCKAAMEHSRLNRIIATV"
elif fragment == "F3" :
    seq = "IHNEVDCDVFFPIDFRSSQSCLPWRKQDHSVLEAWVGSKVPQGKINENGFIYEFEMWIRD"
elif fragment == "total" :
    seq ="DWQKSLTLIVALTLSRGIGLKNDLPWKLKSDMMFFSRVTSGLLVTRSTGQMNVVLMGRKTWESLPAHSRPLKNRINVVISRQEVLDLGGGAYHARSLDDALALLSQIYDSTSKIQLNRVFVIGGGELCKAAMEHSRLNRIIATVIHNEVDCDVFFPIDFRSSQSCLPWRKQDHSVLEAWVGSKVPQGKINENGFIYEFEMWIRD"
else:
    pass

# Set parameters for the rest of the function and generate framework
hm_form_aa = pd.DataFrame(index=aa_list)
hm_form_aa = hm_form_aa.reindex(columns =  list(range(lower,upper)))
mean_list = list()
freq = 0

#Call nature flags
longform["nature"] = ""
for index in list(range(len(longform))):
    if longform.at[index, "codon_WT"] == longform.at[index, "codon_mut"]:
        longform.at[index, "nature"] = "WT"
    elif longform.at[index, "codon_mut"] == "TAG":
        longform.at[index, "nature"] = "stop"
    elif longform.at[index, "aa_WT"] == longform.at[index, "aa_mut"]:
        longform.at[index, "nature"] = "silent"
    else: 
        longform.at[index, "nature"] = "substitution"
        
if condition == "DMSO":
    cols_list = ["T0-1", "T0-2", "T0-3", "DMSO-1", "DMSO-2", "DMSO-3"]
elif condition == "MTX4":
    cols_list = ["T0-1", "T0-2", "T0-3", "MTX4-1", "MTX4-2", "MTX4-3"]
elif condition == "MTX20":
    cols_list = ["T0-1", "T0-2", "T0-3", "MTX20-1", "MTX20-2", "MTX20-3"]


# Correct wild type values with mean of silent mutations
# Here the value that is selected is mutations frequency

silent_df = pd.DataFrame()
for cols in cols_list:
    silent_df[cols] = longform[cols+"_freq"]
silent_df["nature"] = longform["nature"]
silent_df["fragment"] = longform["fragment"]
silent_df["position"] = longform["position"]
silent_df = silent_df[silent_df.nature == "silent"]
silent_df_grouped = silent_df.groupby(by="position").mean()
silent_df_grouped.index = silent_df_grouped.index.astype("int64")
silent_df_grouped["aa"] = ""
list_seq = list(seq)
for posi in silent_df_grouped.index:
    silent_df_grouped.at[posi,"aa"] = list_seq[posi-2]

# Create dataframe to store silent median values
med_sil_df = pd.DataFrame(index=cols_list)
med_sil_df = med_sil_df.reindex(columns = ["F1", "F2", "F3"])


# Calculate median for all samples
for samp in cols_list:
    for fragm in ["F1", "F2", "F3"]:
        sub_df = silent_df[silent_df.fragment == fragm]
        med_samp_frag = np.nanpercentile(sub_df.loc[:,samp], 50)
        med_sil_df.at[samp,fragm] = med_samp_frag
    

# Calculate Selection coefficient for all rows (This is going to be a bitch and a half)
def CoS(longform):
    for ind in list(range(0,len(longform))):
        for rep in [1,2,3]:
            T0 = "T0-"+str(rep)
            T0_freq = "T0-"+str(rep)+"_freq"
            sel = condition+"-"+str(rep)
            sel_freq = condition+"-"+str(rep)+"_freq"
            
            ### Equation to calculate Selection coefficient
            longform.at[ind, sel+"_CoS"] = float((math.log( longform.at[ind,sel_freq] / med_sil_df.at[sel, longform.at[ind,"fragment"] ]) - math.log( longform.at[ind,T0_freq] / med_sil_df.at[T0, longform.at[ind,"fragment"] ]))/longform.at[ind,sel+"_gen"])

CoS(longform)


# Calculate median of CoS
for ind in list(range(len(longform))):
    longform.at[ind, "median_CoS"] = np.median([longform.at[ind, condition+"-1_CoS"], longform.at[ind, condition+"-2_CoS"], longform.at[ind, condition+"-3_CoS"]])

# Save data as final version for selection coefficient for codons
longform.to_csv("./total_"+condition+"_longform_CoS.csv")

#%%
# Visualise data and set T0 read count thresholds using loess regression

# Set condition
condition = "MTX20"

# Load longform data
longform_D = pd.read_csv("./total_DMSO_longform_CoS.csv", index_col=0)
longform_M4 = pd.read_csv("./total_MTX4_longform_CoS.csv", index_col=0)
longform_M20 = pd.read_csv("./total_MTX20_longform_CoS.csv", index_col=0)

# Set proper data based on condition
if condition == "DMSO":
    full_dataframe = longform_D
elif condition == "MTX4":
    full_dataframe = longform_M4
elif condition == "MTX20":
    full_dataframe = longform_M20


# Check different thresholds to identify when Loess regresstion becomes most flat
# This allows to minimize the ffect of low read count on measured selection coeffient.
full_dataframe = full_dataframe[full_dataframe["WT"] != True ]
longform_med = full_dataframe.iloc[:,3:6]
median = lambda x: np.median(x)
medians = longform_med.apply(median, axis = 1)
longform_med["readCount_median"] = medians
longform_med["readCount_median_log2"] = np.log2(longform_med["readCount_median"])
longform_med["median_CoS"] = full_dataframe["median_CoS"]
longform_med["pass_thresh"] = ""
longform_med = longform_med.reset_index(drop = True)

# Selected threshold is 15
thresh = 15

#Add thresholds
for ind in list(range(len(longform_med))):
    medi = longform_med.iloc[ind,0:3]
    if all(i >= thresh for i in medi) == True:
        longform_med.at[ind,"pass_thresh"] = True
    else:
        longform_med.at[ind,"pass_thresh"] = False 

longform_med_all = longform_med.copy()

longform_med = longform_med[["readCount_median_log2","median_CoS", "pass_thresh"]]
longform_med = longform_med[longform_med["pass_thresh"] == True]

sns.lmplot(x="readCount_median_log2", y="median_CoS", data=longform_med,
           lowess=True,scatter_kws={"color": "black"}, line_kws={"color": "C1"});
plt.title("MedianCoS Vs Readcount at T0 in "+condition+" for all - tresh"+str(thresh))
#plt.xlim(0,100)
plt.axvline(np.log2(10), color="black")

plt.savefig('./Figures/loess/CoS_vs_log2medRC_'+condition+'_t'+str(thresh)+'_wt.png', dpi=200, bbox_inches="tight")

#%%
# Generate heatmaps and files based on threshold value

# Generate heatmap aa based on threshold
# =============================================================================
thresh = 15
condition = ""
# =============================================================================
seq ="DWQKSLTLIVALTLSRGIGLKNDLPWKLKSDMMFFSRVTSGLLVTRSTGQMNVVLMGRKTWESLPAHSRPLKNRINVVISRQEVLDLGGGAYHARSLDDALALLSQIYDSTSKIQLNRVFVIGGGELCKAAMEHSRLNRIIATVIHNEVDCDVFFPIDFRSSQSCLPWRKQDHSVLEAWVGSKVPQGKINENGFIYEFEMWIRD"
aa_list = ["G","A","L","M","F","W","K","Q","E","S",
           "P","V","I","C","Y","H","R","N","D","T","*"]

aa_list_prop = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", 
                "A", "V", "I", "L", "M", "F", "Y", "W", "*"]

lower = 2
upper = 206
fragment = "total"


longform_D = pd.read_csv("./total_DMSO_longform_CoS.csv", index_col=0)
#longform_D = longform_D[longform_D.fragment == fragment]
longform_M4 = pd.read_csv("./total_MTX4_longform_CoS.csv", index_col=0)
#longform_M4 = longform_M4[longform_M4.fragment == fragment]
longform_M20 = pd.read_csv("./total_MTX20_longform_CoS.csv", index_col=0)

if condition == "DMSO":
    full_dataframe = longform_D
elif condition == "MTX4":
    full_dataframe = longform_M4
elif condition == "MTX20":
    full_dataframe = longform_M20


# Runs fresh analysis from previous cells
longform_med = full_dataframe.iloc[:,3:6]
median = lambda x: np.median(x)
medians = longform_med.apply(median, axis = 1)
longform_med["readCount_median"] = medians
longform_med["readCount_median_log2"] = np.log2(longform_med["readCount_median"])
longform_med["median_CoS"] = full_dataframe["median_CoS"]
longform_med["pass_thresh"] = ""
longform_med = longform_med.reset_index(drop = True)


#Add thresholds
for ind in list(range(len(longform_med))):
    medi = longform_med.iloc[ind,0:3]
    if all(i >= thresh for i in medi) == True:
        longform_med.at[ind,"pass_thresh"] = True
    else:
        longform_med.at[ind,"pass_thresh"] = False 

passed = longform_med["pass_thresh"]
full_dataframe["pass_thresh"] = passed
longform = full_dataframe
longform = longform[longform["pass_thresh"] == True]
longform = longform.reset_index(drop = True)

hm_form_aa_CoS = pd.DataFrame(index=aa_list_prop)


# Aggregate all aa per codon with selection coefficients using codon median
hm_form_aa_CoS = hm_form_aa_CoS.reindex(columns =  list(range(lower,upper)))
mean_list = list()
freq = 0

for pos in list(range(lower,upper)):
    for aa in aa_list_prop:
        for line in list(range(0,len(longform.index))):
            if longform.at[line,'position'] == pos and longform.at[line,'aa_mut'] == aa:
                mean_list.append(float(longform.at[line,'median_CoS']))
            else:
                pass
        if len(mean_list) == 0:
            pass
        else:
            freq = statistics.median(mean_list)
            #print(aa+"_"+str(freq))
            mean_list = list()
            hm_form_aa_CoS.at[aa,pos] = freq
            freq = 0    

# Run function on dataframe to create heatmap format
#Create "safety" copy of DF
hm_form_aa_pre_wt = hm_form_aa_CoS.copy()



# Get meadian CoS for silent mutations to fill WT
silent_df_CoS = pd.DataFrame()
silent_df_CoS["median_CoS"] = longform["median_CoS"]
silent_df_CoS["fragment"] = longform["fragment"]
silent_df_CoS["nature"] = longform["nature"]
silent_df_CoS["position"] = longform["position"]
silent_df_CoS = silent_df_CoS[silent_df_CoS.nature == "silent"]
silent_df_grouped_CoS = silent_df_CoS.groupby(by="position").mean()
silent_df_grouped_CoS.index = silent_df_grouped_CoS.index.astype("int64")


silent_df_grouped_CoS["aa"] = ""
list_seq = list(seq)
for posi in silent_df_grouped_CoS.index:
    silent_df_grouped_CoS.at[posi,"aa"] = list_seq[posi-2]

for posit in silent_df_grouped_CoS.index:
    hm_form_aa_CoS.at[silent_df_grouped_CoS.at[posit,"aa"], posit] = 0

index = lower
for codon in list(seq):
    value = hm_form_aa_CoS.at[codon, index]
    isNan = np.isnan(value)
    if isNan == True:
        if index <= 74:
            hm_form_aa_CoS.at[codon, index] = 0
            index = index + 1
        elif index > 74 and index <= 146:
            hm_form_aa_CoS.at[codon, index] = 0
            index = index + 1
        else:
            hm_form_aa_CoS.at[codon, index] = 0
            index = index + 1
    else:
        index = index + 1
        

#longform.to_csv("./Threshold/"+fragment+"_"+condition+"_longform_CoS_t"+str(thresh)+".csv")
hm_form_aa_CoS.to_csv("./Threshold/"+fragment+"_"+condition+"_heatmap_format_CoS_t"+str(thresh)+"_aa-prop.csv")

# Dataframe with all selection coefficients are read for downstream analysis. 
