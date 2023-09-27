#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 10:26:39 2023

@author: francois
"""
import os
import subprocess
from Bio.Seq import Seq
#%%
# Generate all possible single aa mutants in nt from nt sequence
# Must be a multiple of 3, since permutations are done on a codon basis using NNN
# Output is fasta format file with all possible mutant nt sequences
### Header format is: wt codon ; wt aa ; position ; mutant aa ; mutant codon
# =============================================================================
 # Sequence in nt
seq = "ATGGACTGGCAAAAGTCCTTAACTTTGATCGTTGCTTTGACTTTGTCTAGAGGTATTGGTTTGAAAAACGACTTGCCTTGGAAGTTGAAGTCTGACATGATGTTCTTTTCCAGAGTTACTTCTGGTCTGTTGGTCACCAGATCAACCGGTCAAATGAACGTCGTTTTGATGGGTCGTAAGACCTGGGAATCTTTACCAGCTCATTCTAGACCATTAAAGAACAGAATCAACGTTGTCATCTCCAGACAAGAAGTCTTGGATTTAGGTGGTGGTGCTTACCACGCTAGATCTTTGGACGATGCCCTAGCCTTGTTGTCCCAAATCTACGATTCTACTTCTAAAATCCAATTGAACAGAGTTTTCGTTATCGGTGGCGGTGAATTGTGTAAGGCTGCTATGGAACACAGTAGATTGAACCGTATCATTGCCACCGTCATTCACAACGAAGTTGACTGTGATGTCTTCTTCCCAATTGACTTCAGATCTTCTCAATCCTGTTTGCCATGGAGAAAGCAAGATCACTCCGTCTTGGAAGCTTGGGTTGGTTCCAAGGTTCCACAAGGTAAGATTAATGAAAACGGTTTCATTTACGAATTCGAAATGTGGATTCGTGACATCtaa"
 # What position if the first codon of your sequence in your protein (only use if using fragments. If not, default is 1)
index_pos = 1
 # Output file name that you want
output_name = "F3_all_variants.fasta"
# =============================================================================

# Get length in codons
nt_len = len(seq)
# Get length in aa
aa_len = int(len(seq)/3)
# Split sequence in "codons"
list_seq_codon = [seq[i:i+3] for i in range(0, len(seq), 3)]

# Create dictionary of nt to aa
gencode = { 'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
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
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W' }

#Create list of codons
codons = list(gencode.keys())
# Index a list of the length of your sequence
positions = list(range(0,aa_len))

pos = 199
print(list_seq_codon[pos-2],list_seq_codon[pos-1],list_seq_codon[pos])


# Create all possible sequences based on the input sequence
## Open file to write in 
with open(output_name, 'w') as source:
    # For all positions in the sequence
    for pos in positions:
        # Find WT codon in your sequence
        wt_codon = list_seq_codon[pos]
        # Find wt aa
        wt_aa = gencode[wt_codon]
        # For all possible subsitutitons
        for codon in codons:
            # Set identitiy of new codon
            mut_codon = codon
            # Set meta information (Wt codon; Wt aa; positions; mutant aa; mutant codon)
            mutseq_id = ">"+wt_codon+";"+str(wt_aa)+";"+str(pos+index_pos)+";"+gencode[codon]+";"+codon+"\n"
            list_mutseq = list()
            list_mutseq = list_seq_codon[:pos]
            list_mutseq.append(mut_codon)
            list_mutseq.extend(list_seq_codon[pos+1:])
            mutseq = "".join(list_mutseq)
            source.write(mutseq_id)
            source.write(mutseq+"\n")
        
source.close()


#%%
# As function

def all_mutants(nt_seq, output_name, index_pos = 1):
    # Get length in codons
    nt_len = len(nt_seq)
    # Get length in aa
    aa_len = int(len(nt_seq)/3)
    # Split sequence in "codons"
    list_seq_codon = [nt_seq[i:i+3] for i in range(0, len(nt_seq), 3)]
    
    # Create dictionary of nt to aa
    gencode = { 'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
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
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W' }
    
    #Create list of codons
    codons = list(gencode.keys())
    # Index a list of the length of your sequence
    positions = list(range(0,aa_len))
    
    
    # Create all possible sequences based on the input sequence
    ## Open file to write in 
    with open(output_name, 'w') as source:
        # For all positions in the sequence
        for pos in positions:
            # Find WT codon in your sequence
            wt_codon = list_seq_codon[pos]
            # Find wt aa
            wt_aa = gencode[wt_codon]
            # For all possible subsitutitons
            for codon in codons:
                # Set identitiy of new codon
                mut_codon = codon
                # Set meta information (Wt codon; Wt aa; positions; mutant aa; mutant codon)
                mutseq_id = ">"+wt_codon+";"+str(wt_aa)+";"+str(pos+index_pos)+";"+gencode[codon]+";"+codon+"\n"
                list_mutseq = list()
                list_mutseq = list_seq_codon[:pos]
                list_mutseq.append(mut_codon)
                list_mutseq.extend(list_seq_codon[pos+1:])
                mutseq = "".join(list_mutseq)
                source.write(mutseq_id)
                source.write(mutseq+"\n")
            
    source.close()

