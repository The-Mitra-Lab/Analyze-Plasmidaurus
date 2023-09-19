# Analyze-Plasmidsaurus

<p align="center">
    <img src="https://github.com/The-Mitra-Lab/Analyze-Plasmidaurus/blob/main/analyze_plasmidsaurus2.png">
</p>

**Analyze-plasmidsaurus** is a python shell script developed by the [Mitra Lab](https://mitralab.wustl.edu/) to annotate mutations in plasmidsaurus GenBank files. 

This script takes a reference sequence in file GenBank format and plasmidsaurus file  in GenBank format and conducts a global alignment between the two sequences. Differences between the two sequences are tracked and the type of sequence mutation is added as an annotation to the plasmidsaurus Genbank file. Base pair changes and subsequent Amino Acid changes are included in the annotation.    \

Once annotations are added to plasmidsaurus GenBank file count_mutations.py can be used to count frameshifts, substitutions, transitions, transversions, and coding changes. The counts are outputed as csv file. 


## Mutation Annotations

**Substitution; AA change:** A->T; Pro -> Val \
**Substitution; Silent mutaiton:** A->T; Silent \
**Frameshift; Deletion:** A-> - ; Pro -> Deletion \
**Frameshift; Insertion:** - -> T; Insertion -> Val 

## Script Usage 

**To annotate a single GenBank file use:** \
python3 analyze_plasmidsaurus.py reference.gbk plasmidsaurus.gbk library_name

**To annotate multiple GenBank files use:** \
for FILE in directory\*.gbk; do python3 analyze_plasmidsaurus.py directory\reference.gbk $FILE library_name;Done 

**To count mutations in annotated GenBank file use:** \
python3 count_mutations.py annotated_file.gbk clone_name library_name 

**To count mutations in annotated multiple GenBank (in same library) file use:**    \
for FILE in directory\*annotated.gbk; do python3 count count_mutations.py $FILE;Done 




