#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 14:04:51 2023

@author: francois
"""
# PyMol scripts

# These functions are used to calculate the distance matrix between
# all residues of two chains. This is useful to calculate distance 
# between ligand and all residues of a protein/chain.



# Import necessary packages and libraries
from Bio import PDB
import Bio.PDB
import numpy
import pandas as pd
import os
aa_list_nostop = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", 
                "A", "V", "I", "L", "M", "F", "Y", "W"]
#%%
# Define your protein's parameters and create the python object
pdb_code = "Pj_DHFR"
pdb_filename = "/Users/francois/Dropbox/Mac/Desktop/Pj_pred_MTX_NDP_dist.pdb"
structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)
model = structure[0]


#%%
# WT function - Gives distance from all CA of chain one to all CA of chain two

## Here, define your chain one as being you ligand. Make sure that you chose 
## a CA for this ligand as not all ligand will have a CA

def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return numpy.sqrt(numpy.sum(diff_vector * diff_vector))

def calc_dist_matrix(chain_one, chain_two) :
    """Returns a matrix of C-alpha distances between two chains"""
    answer = numpy.zeros((len(chain_one), len(chain_two)), numpy.float)
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            answer[row, col] = calc_residue_dist(residue_one, residue_two)
    return answer

dist_matrix = calc_dist_matrix(model["X"], model["A"])


#%%
# Min distance function- Gives minimal distance between ligand and residue CA

## Similar to above, but calculates distance from all atoms of a residue of 
## chain/molecule one to all CA of molecule two. This is useful when chain one is
## larger, and/or not a protein. I used this for my different ligands

def calc_residue_dist_min(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    atom_min = pd.Series(dtype="float64")
    for ro, atom in enumerate(residue_one):
        atom = str(atom)
        inte = atom.split(" ")
        atom = inte[1].split(">")[0]
        if atom[0].startswith("H"):
            pass
        else:
            diff_vector  = residue_one[atom].coord - residue_two["CA"].coord
            atom_min.loc[ro] = numpy.sqrt(numpy.sum(diff_vector * diff_vector))
            min_dist = atom_min.min()
    return min_dist

def calc_dist_matrix_min(chain_one, chain_two) :
    """Returns a matrix of C-alpha distances between two chains"""
    answer = numpy.zeros((len(chain_one), len(chain_two)), float)
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            answer[row, col] = calc_residue_dist_min(residue_one, residue_two)
    return answer


dist_matrix_min_CA = calc_dist_matrix_min(model["X"], model["A"])

#%%
# Plotting the distance map
import pylab
pylab.matshow(numpy.transpose(dist_matrix))
pylab.colorbar()
pylab.show()

#%%
dist_matrix_min = answer.copy()
dist_matrix_min_t = dist_matrix_min.transpose()

for re, val in enumerate(dist_matrix_min_t):
    print(val[0])


#%%

# ALL VS ALL


runcell(0)
runcell(1)
tot_dist_df = pd.DataFrame(columns = ("pos_residue", 'pos_ligand', 'atom_ligand', 'index_resdue', 'atom_resdue', 'distance'))

## DEBUG
chain_one = model["X"]
chain_two = model["A"]
ligand_atoms = pd.Series(dtype="float64")

#"""Returns tall atoms of chain vs all atoms of ligand distance"""
for row, residue_one in enumerate(chain_one) :
    for col, residue_two in enumerate(chain_two) :
        atom_min = pd.Series(dtype="float64")
        res_min = pd.Series(dtype="float64")
        for ro, atom in enumerate(residue_one):
            inte = str(atom).split(" ")
            atom = inte[1].split(">")[0]
            if atom[0].startswith("H"):
                pass
            else:
                for ro2, atom2 in enumerate(residue_two):
                    atom2 = str(atom2)
                    inte = atom2.split(" ")
                    atom2 = inte[1].split(">")[0]
                    if atom2[0].startswith("H"):
                        pass
                    elif len(atom2) == 1:
                        pass
                    else:
                        #print(col, ro, atom, ro2, atom2)
                        diff_vector_res = residue_one[atom].coord - residue_two[atom2].coord
                        info = [col, ro, atom, ro2, atom2, numpy.sqrt(numpy.sum(diff_vector_res * diff_vector_res))]
                        tot_dist_df.loc[len(tot_dist_df)] = info

#%%

conv = {'distance': float, 'pos_residue': int, 'pos_ligand': float}
tot_dist_df = tot_dist_df.astype(conv)
print(tot_dist_df.dtypes)

tot_dist_df_test = tot_dist_df.copy()
tot_dist_df_test['pos_residue'] = tot_dist_df_test['pos_residue'] +1
#tot_dist_df_test = tot_dist_df_test.groupby('pos_residue', as_index=False, group_keys=False).min()



min_distances = tot_dist_df_test.groupby('pos_residue')['distance'].min()

# Create a new dataframe using the minimum distances
new_df = tot_dist_df_test.loc[tot_dist_df_test.groupby('pos_residue')['distance'].idxmin()]
new_df.to_csv('./WT_mindist_allatoms.csv')
#%%
#"""Returns tall atoms of chain vs all atoms of ligand distance"""
# ALL IN FLEXDDG
runcell(0)

#mutation = 'Y'
#tot_dist_df = pd.DataFrame(columns = ('nstruct','residue', "pos_residue", 'pos_ligand', 'atom_ligand', 'index_resdue', 'atom_resdue', 'distance'))
#rootdir = '/Users/francois/Dropbox/Mac/Desktop/strucutre_colormap/FlexddG_199/MTX/output_saturation100/Pj_DFR_A199_'+mutation+'/'

for mutation in aa_list_nostop:
    tot_dist_df = pd.DataFrame(columns = ('nstruct','residue', "pos_residue", 'pos_ligand', 'atom_ligand', 'index_resdue', 'atom_resdue', 'distance'))
    rootdir = '/Users/francois/Dropbox/Mac/Desktop/strucutre_colormap/FlexddG_199/MTX/output_saturation100/Pj_DFR_A199_'+mutation+'/'
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            cur_dir = os.path.join(subdir, file)
            information = cur_dir.split('/')
            if len(information) < 13:
                pass
            else:
                if information[12] == 'mut_35000.pdb':
                    pdb_code = "Pj_DHFR"
                    structure = Bio.PDB.PDBParser().get_structure(pdb_code, cur_dir)
                    model = structure[0]
                    chain_one = model["X"]
                    chain_two = model["A"]
                    for row, residue_one in enumerate(chain_one) :
                        ligand = residue_one.__repr__().split(' ')
                        if ligand[1] == 'MTX':
                            for col, residue_two in enumerate(chain_two) :
                                if col != 198:
                                    pass
                                else:
                                    for ro, atom in enumerate(residue_one):
                                        inte = str(atom).split(" ")
                                        atom = inte[1].split(">")[0]
                                        if atom[0].startswith("H"):
                                            pass
                                        else:
                                            for ro2, atom2 in enumerate(residue_two):
                                                atom2 = str(atom2)
                                                inte = atom2.split(" ")
                                                atom2 = inte[1].split(">")[0]
                                                if atom2[0].startswith("H"):
                                                    pass
    # =============================================================================
                                                elif len(atom2) == 1:
                                                    pass
    # =============================================================================
                                                elif atom2[0].isdigit(): 
                                                    pass
                                                else:
                                                    diff_vector_res = residue_one[atom].coord - residue_two[atom2].coord
                                                    aa = residue_two.__repr__().split(' ')[1]
                                                    info = [information[11], aa, col+1, ro, atom, ro2, atom2, numpy.sqrt(numpy.sum(diff_vector_res * diff_vector_res))]
                                                    tot_dist_df.loc[len(tot_dist_df)] = info
                    
    
    
    
    conv = {'distance': float, 'pos_residue': int, 'pos_ligand': float}
    tot_dist_df = tot_dist_df.astype(conv)
    
    tot_dist_df_test = tot_dist_df.copy()
    tot_dist_df_test['pos_residue'] = tot_dist_df_test['pos_residue']
    #tot_dist_df_test = tot_dist_df_test.groupby('pos_residue', as_index=False, group_keys=False).min()
    
    
    
    min_distances = tot_dist_df_test.groupby('pos_residue')['distance'].max()
    
    # Create a new dataframe using the minimum distances
    new_df = tot_dist_df_test.loc[tot_dist_df_test.groupby('nstruct')['distance'].idxmin()]
    new_df.to_csv('./'+mutation+'_mindist_allatoms.csv')
    
    
    
    











