#!/bin/sh

#usage: analyze_plasmidsaurus.py "reference_genome.gbk" "plasmidsaurus.gbk" "library_name"



import glob
import os
import sys



from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import GenBank
from Bio import SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation, FeatureLocation
from itertools import groupby
from operator import itemgetter




# tf name needs to be *TF_*
ref_genome = sys.argv[1]
plasmid_genome = sys.argv[2]
lib_name = sys.argv[3]
file_name = os.path.splitext(plasmid_genome)




def translate(seq):
    """
          Translates nucelotide sequence to aa sequence and returns translated aa sequence and start and end positions of aa

    """
    codontable = {'ATA': 'Ile', 'ATC': 'Ile', 'ATT': 'Ile', 'ATG': 'Met',
              'ACA': 'Thr', 'ACC': 'Thr', 'ACG': 'Thr', 'ACT': 'Thr',
              'AAC': 'Asn', 'AAT': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
              'AGC': 'Ser', 'AGT': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
              'CTA': 'Leu', 'CTC': 'Leu', 'CTG': 'Leu', 'CTT': 'Leu',
              'CCA': 'Pro', 'CCC': 'Pro', 'CCG': 'Pro', 'CCT': 'Pro',
              'CAC': 'His', 'CAT': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
              'CGA': 'Arg', 'CGC': 'Arg', 'CGG': 'Arg', 'CGT': 'Arg',
              'GTA': 'Val', 'GTC': 'Val', 'GTG': 'Val', 'GTT': 'Val',
              'GCA': 'Ala', 'GCC': 'Ala', 'GCG': 'Ala', 'GCT': 'Ala',
              'GAC': 'Asp', 'GAT': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
              'GGA': 'Gly', 'GGC': 'Gly', 'GGG': 'Gly', 'GGT': 'Gly',
              'TCA': 'Ser', 'TCC': 'Ser', 'TCG': 'Ser', 'TCT': 'Ser',
              'TTC': 'Phe', 'TTT': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu',
              'TAC': 'Tyr', 'TAT': 'Tyr', 'TAA': '*', 'TAG': '*',
              'TGC': 'Cys', 'TGT': 'Cys', 'TGA': '*', 'TGG': 'Trp',
              '---': '-',
              }
    aa_list = []
    coord_start = []
    coord_end = []

    for i in range(0, len(seq), 3):
        if seq[i:i + 3] in codontable:
            codon = seq[i:i + 3]
            aa_list.append(codontable[codon])
            coord_start.append(i)
            coord_end.append(i+2)
#
        else:
            residue = "X"
            aa_list.append(residue)
            coord_start.append(i)
            coord_end.append(i+2)


#

    return coord_start, coord_end, aa_list


def get_aa(start, end, aa):
    """
    Returns list of amino acids and a list of amino acid cooridinate where a subsitution has occured in nucleotide sequence
    """
    aa_list = []
    aa_coord = []
    for x,y,z in zip(start, end, aa):
        for i in filtered_coord:
            if i in range(x,y+1):
                aa_list.append(z)
                aa_coord.append(i)
    return aa_list, aa_coord

def add_feature(start, end, plasmid_record, hex_code):
    """
    can add subsitution, deletion, or insertion feature to genbank file
    """
    # Use the locations do define a FeatureLocation
    my_feature_location2= FeatureLocation(start, end)
    # Define a feature type as a text string
    my_feature_type2 = "misc_feature"
    # Create a SeqFeature
    my_feature2 = SeqFeature(my_feature_location2,type=my_feature_type2,qualifiers={'label':[f], 'note':['pLannotate', hex_code]})
    #  Append your newly created SeqFeature to your SeqRecord
    plasmid_record.features.append(my_feature2)


plasmid_g = SeqIO.read(plasmid_genome, "genbank")
ref_g = SeqIO.read(ref_genome, "genbank")


#filteres out plasmid files that are 50 bp longer or shorter than reference plasmid files
#Cutoff can be less (~20bp)
if len(ref_g) - len(plasmid_g) < 50:

#changes the ori of the plasmid sequence to be the same as the reference
    ori_seq = ref_g.seq[0:25]
    ori_coordinate = []
    if ori_seq in plasmid_g:
        ori_coordinate.append(plasmid_g.seq.find(ori_seq))
    else:
        reverse_seq = plasmid_g.seq.reverse_complement()
        ori_coordinate.append(reverse_seq.find(ori_seq))
        plasmid_g.seq = reverse_seq

    shifted = plasmid_g[ori_coordinate[0]:] + plasmid_g[:ori_coordinate[0]]


    new_plasmid_record = SeqRecord(shifted)
    new_plasmid_record.features = ref_g.features
    new_plasmid_record.annotations = ref_g.annotations

#aligns the plasmid sequence to the reference sequence
#saves the aligned sequence as ref_seq and plasmid_seq

    alignments = pairwise2.align.globalms(ref_g.seq, shifted.seq, 5, -4, -4, -.1)
    formated_seq = format_alignment(*alignments[0])
    ref_seq = formated_seq.split("\n")[0]
    plasmid_seq = formated_seq.split("\n")[2]

#if there is a difference between the two sequences it appends the coordinate to a list
    coord = [i for i in range(len(ref_seq)) if ref_seq[i] != plasmid_seq[i]]

# this only keeps the coordinates within the gene of interest
#this can be removed if you are intrested in all mutations in the plasmid sequence
    filtered_coord = []
    for feature in new_plasmid_record.features:
        if feature.type == "gene":
            for x in coord:
                if x > feature.location.start and x < feature.location.end:
                    filtered_coord.append(x)



    lables = []
    for location in filtered_coord:
        lables.append(ref_seq[location]+"->"+plasmid_seq[location])



#translates the reference and plasmid sequences
    ref_start, ref_end, ref_aa  = translate(ref_seq)

    plasmid_start, plasmid_end, plasmid_aa  = translate(plasmid_seq)

    ref_aa_list, ref_coord = get_aa(ref_start, ref_end, ref_aa)



    plasmid_aa_list, plasmid_coord = get_aa(plasmid_start, plasmid_end, plasmid_aa)


    s = []
    e = []
    big_frameshift = []

#finds coordinates that are consecutive

    for k, g in groupby(enumerate(plasmid_coord), lambda ix : ix[0] - ix[1]):

            groups = list(map(itemgetter(1), g))

            if len(groups) > 1:

                big_frameshift.append(groups)
                s.append(groups[0])
                e.append(groups[-1])
            else:
                pass

    filtered_aa_coord = []
    filtered_aa_ref = []
    filtered_aa_plasmid = []
    frame_shift_coord = []
    frame_ref_aa = []
    frame_plasmid_aa = []

#puts coordinates associated with a large frame shift into seperate list
    for a,b,c,d in zip(ref_aa_list, ref_coord, plasmid_aa_list, plasmid_coord):
        if not big_frameshift:
            filtered_aa_coord.append(b)
            filtered_aa_ref.append(a)
            filtered_aa_plasmid.append(c)
        else:
            if len(big_frameshift) > 1:
                if b in big_frameshift and a == "X" or a== "-":
                    frame_shift_coord.append(b)
                    frame_ref_aa.append(a)
                    frame_plasmid_aa.append(c)
                elif d in coord and c == "X" or c == "-":
                    frame_shift_coord.append(d)
                    frame_ref_aa.append(a)
                    frame_plasmid_aa.append(c)
                else:
                    filtered_aa_coord.append(b)
                    filtered_aa_ref.append(a)
                    filtered_aa_plasmid.append(c)
            else:
                for coord in big_frameshift:
                    if b in coord and a == "X" or a== "-":
                        frame_shift_coord.append(b)
                        frame_ref_aa.append(a)
                        frame_plasmid_aa.append(c)
                    elif d in coord and c == "X" or c == "-":
                        frame_shift_coord.append(d)
                        frame_ref_aa.append(a)
                        frame_plasmid_aa.append(c)
                    else:
                        filtered_aa_coord.append(b)
                        filtered_aa_ref.append(a)
                        filtered_aa_plasmid.append(c)

#replaced 'X' with either 'Insertion' or 'Deletion'
    filtered_aa_ref = list(map(lambda filtered_aa_ref: filtered_aa_ref.replace('X', 'Insertion'), filtered_aa_ref))
    filtered_aa_plasmid = list(map(lambda filtered_aa_plasmid: filtered_aa_plasmid.replace('X', 'Deletion'), filtered_aa_plasmid))


#creates lables for aa change
    lables = []
    for location in filtered_aa_coord:
        lables.append(ref_seq[location]+"->"+plasmid_seq[location])
    aa_change = []
    for ref, plas in zip(filtered_aa_ref, filtered_aa_plasmid):
        if ref == plas:
            aa_change.append('silent')
        else:
            aa_change.append(ref+"->"+plas)

    ref_list = []
    plasmid_list = []


#creates a lable for large frameshift
    for location in big_frameshift:
        ref_string = ""
        plasmid_string = ""
        for loc in location:
            # print(loc)
            if ref_seq[loc] == '-':
                ref_string += '-'
            else:
                ref_string += ref_seq[loc]
            if plasmid_seq[loc] == '-':
                plasmid_string  += '-'
            else:
                plasmid_string += plasmid_seq[loc]
        ref_list.append(ref_string)
        plasmid_list.append(plasmid_string)

    frameshift_lable = []
    for ref, plasmid in zip(ref_list, plasmid_list):
        if '-' in ref and '-' not in plasmid:
            frameshift_lable.append('Insertion'+ '->'+ plasmid)
        else:
            frameshift_lable.append(ref + '->' + 'Deletion')

    zp = list(zip(lables, aa_change))

#creates a new seq record and adds annotations from reference record
    new_plasmid_record = SeqRecord(shifted.seq)
    new_plasmid_record.features = ref_g.features
    new_plasmid_record.annotations = ref_g.annotations

#adds subsitution features and single deletions and insertions
    for f, b in zip(zp, filtered_aa_coord):
        if 'Deletion' in f[1]:
            add_feature(b, b, new_plasmid_record, 'color: #ee10ec')
            
        elif 'Insertion' in f[1]:

            add_feature(b, b+1, new_plasmid_record, 'color: #ecee10')

        else:

            add_feature(b, b+1, new_plasmid_record, 'color: #ff0000')

#adds features for large frameshift
    for start_pos, end_pos, lable in zip(s, e, frameshift_lable):
        if 'Deletion' in lable:
            add_feature(start_pos, start_pos, new_plasmid_record, 'color: #ee10ec')
        else:
            add_feature(start_pos, end_pos, new_plasmid_record, 'color: #ecee10')



    file_name = file_name[0].replace('pLann',"")
    output_file = open(file_name+lib_name + '.gbk', 'w')
    SeqIO.write(new_plasmid_record, output_file, 'genbank')
