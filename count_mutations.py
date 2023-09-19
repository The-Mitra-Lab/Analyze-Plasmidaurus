#!/bin/sh

#usage: python3 count_mutations.py "file.gbk" "lib_name"
#This scipt will output a csv file that with frameshift, subs, transitions, treansversions, and coding change counts 

import glob
import os
import sys
import pandas as pd


from Bio import SeqIO
from Bio import GenBank
from Bio import SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SimpleLocation




from snapgene_reader import snapgene_file_to_dict, snapgene_file_to_seqrecord




filepath = sys.argv[1]
lib_name = sys.argv[2]


record = SeqIO.read(filepath, "genbank")

f,t  = os.path.splitext(os.path.basename(filepath))



mut = []
for seqrecord in record.features:
    for subList in seqrecord.qualifiers.values():
        for s in subList:
            if '->' in s:
                mut.append(s)


mut =  mut[2:]


transition_count = 0
transversion_count = 0
frame_shift = 0
for mutation in mut:
    if 'A->T'in mutation:
        transition_count+=1
    elif 'T->A'in mutation:
        transition_count+=1
    elif 'G->C' in mutation:
        transition_count+=1
    elif 'C->G'in mutation:
        transition_count+=1
    elif 'Insertion' in mutation:
        frame_shift+=1
    elif 'Deletion' in mutation:
        frame_shift+=1
    else:
        transversion_count+=1



sub_count = transition_count + transversion_count

silent_mut = 0
coding_change = 0
for mut_type in mut:
    if "Deletion" in mut_type:
        pass
    elif "Insertion" in mut_type:
        pass
    elif mut_type[10] == mut_type[13]:
        silent_mut+=1
    else:
        coding_change+=1


data = {'Clone':[f], 'Frameshift':[frame_shift], 'Subs':[sub_count], 'Transition':[transition_count], 'Transversion':[transversion_count], 'Coding Change':[coding_change]}
df = pd.DataFrame(data)

with open(lib_name+ "_summary.csv","a") as f:
    df.to_csv(f,header=False,index=False)
