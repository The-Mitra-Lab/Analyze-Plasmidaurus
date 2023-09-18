## Analyze-Plasmidsaurus

<p align="center">
    <img src="https://github.com/The-Mitra-Lab/Analyze-Plasmidaurus/blob/main/analyze_plasmidsaurus.png">
</p>

**Analyze-plasmidsaurus** is a python script developed by the [Mitra Lab](https://mitralab.wustl.edu/) to annotate mutations in plasmidsaurus GenBank files. 

This script takes a reference sequence in file GenBank format and plasmidsaurus file  in GenBank format and conducts a global alignment between the two sequences. Differences between the two sequences are tracked and the type of sequence mutation is added as an annotation to the plasmidsaurus Genbank file. Basepair changes and subsequent Amino Acid changes are included in the annotation. 


# Mutation Annotations

**Substitution; AA change:** A->T; Proline -> Valine \
**Substitution; Silent mutaiton:** A->T; Silent \
**Frameshift; Deletion:** A-> - ; Proline -> Deletion \ 
**Frameshift; Insertion:** - -> T; Insertion -> Valine \  

