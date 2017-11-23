# RothSchulze_microbiome

### Current pipeline use for 16s rRNA sequence processing for the microbiome analysis in the Environmental Determinant of Islet Autoimmunity (ENDIA) project

The script "Microbiome 16S rRNA amplicon sequences processing pipeline.sh" contains a set instructions on programs and their commands that can be used to process raw 16S rRNA sequences to obtain a biom file that contains the OTU table, taxonomic information and metadata. Additionally some instructions are included to build a phylogenetic tree with the representative sequences from each OTU and also to normalize the OTU table with the Cumulative Sum Scaling (CSS) method.

Requirements to run this pipeline are to install:

QIIME versions 1.8.0 and 1.9.1
PEAR v0.9.0 
usearch version 9.2
vsearch v2.4.3 
MOTHUR v.1.38.1
perl v5.14.0

The resulting biom file can be used as input for other microbiome analysis programs such as phyloseq.



