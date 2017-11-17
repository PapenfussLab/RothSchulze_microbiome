## Microbiome 16S rRNA amplicon sequences processing pipeline 
## 2016-2017 Alexandra J. Roth Schulze and Jocelyn Sietsma Penington

## Used: QIIME (versions 1.8.0 and 1.9.1), PEAR v0.9.0, usearch (version 9.2), 
# vsearch v2.4.3, MOTHUR v.1.38.1. and perl v5.14.0

## Note: $PATH denotes the directory in which your files are contained or where you would
## like to have your output. Therefore, you should change this at your covenience

## 1.- Quality inspection of raw data:
fastqc -o $PATH/FolderContainingFastqFiles --noextract $PATH/*fastq.gz 

# Uncompress
gzip -cd $PATH/MISEQXXXX_Pair1.fastq.gz > $PATH/MISEQXXXX_R1.fastq
gzip -cd $PATH/MISEQXXXX_Pair2.fastq.gz > $PATH/MISEQXXXX_R2.fastq

## 2.- Creating and checking the mapping file:

# You need to start with your sample metadata mapping file which contains the per-sample 
# barcode sequences and other technical information. For more info and examples visit: 
# http://qiime.org/documentation/file_formats.html

# Check if the mapping file has errors
validate_mapping_file.py -m $PATH/Mapping_file.csv

## 3.- Extracting the barcodes:
extract_barcodes.py -f ../1-Raw_data_and_quality/MISEQ2128_S1_L001_R1_001.fastq  \
-r $PATH/MISEQXXXX_R1.fastq -c barcode_paired_end -o $PATH/bar_exed_sep_ends -l 8 -L 8 \
-m $PATH/Mapping_file.csv --attempt_read_reorientation

# Ouput files:
# barcodes.fastq  barcodes_not_oriented.fastq  reads1.fastq  reads2.fastq  
# Note: If the resulting file reads1.fastq has in the header a 1, the order of the primers 
# in the mapping file is correct

## 4.- Merge overlapping paired-end reads with PEAR.
# Options: -v is min overlap, -m is max assembled length, -n is min assembled length
# -j is number of threads
pear -f $PATH/MISEQXXXX_R1.fastq -r $PATH/MISEQXXXX_R2.fastq  -v 50 -m 600 -n 300 \
-j $(nproc) -o $PATH/MISEQXXXX

# Output files:
# $PATH/MISEQXXXX.assembled.fastq MISEQXXXX.assembled.fastq

# Note: PEAR essentially add the 2 quality scores if the calls agree (with mathematical justification). 
# However, this breaks the conventions of the phred score in FASTQ. It doesn't break FASTQ completely
# - values up to ASCII 126 = Phred+33 93 - can be used, but they are unconventional! 
# But this is something QIIME 1.9.0’s split_libraries_fastq.py (for demultiplexing step 5) 
# cannot handle and therefore QIIME 1.8.0 must be used for the split_libraries_fastq.py step

## 4.- Discard sequences in barcodes.fastq that are not in sequences file:
# For this step a python script written by Jocelyn Penington available it GitHub is used 

python $PATH/trim_fastq_to_matching.py -f bar_exed.assembled.fastq \
 -m bar_exed_sep_ends/barcodes.fastq  -o MISEQXXXX.barcodematched.fastq
 
## 5.- Demultiplexing fastq sequencing data:
# split_libraries_fastq: label sequences with sample ID based on index sequences.

# You need to use Qiime 1.8.0 for this step

lamboot

# Options: -p --min_per_read_length_fraction, -q --phred_quality_threshold, -n --sequence_max_n
  
split_libraries_fastq.py  -i $PATH/bar_exed/reads.fastq -b $PATH/bar_exed/barcodes.fastq \
-m $PATH/mapping.txt --barcode_type 16 -q 29 -p 0.90 --phred_offset 33 -n 1 -o $PATH/labelled_hiqual 

# Output: histograms.txt  seqs.fna  split_library_log.txt

## 6.-  Remove universal primers and amplicon primers
# For this a script written by Jocelyn Penington available in GitHub is used. This trims 
# 16S rRNA primers as well as Illumina universal sequencing primers
python $PATH/trim_fasta_amplicons.py -i seqs.fna -o trimmed_seqs.fna

## 7.- Align sequences and cut the alignment (within MOTHUR):
# Note: This is something implemented in MOTHUR that for some reason is not done in QIIME. 
# The reason to do this is to keep only sequences from the same region of the 16S rRNA gene
# and to have all the reads with the exactly same length (meaning exactly the same region!).
# For this we use the silva database from MOTHUR:

mothur
 
# A) Align to silva.bacteria database (you should copy the database to your directory).
align.seqs(fasta=trimmed_seqs.fna, reference=silva.bacteria.fasta, flip=t, processors=24)
# Output: trimmed_seqs.align

# B) Check in which bases are most of the sequences aligned
summary.seqs(fasta=trimmed_seqs.align, processors=24)

# Output Example:
#                 Start   End     NBases  Ambigs  Polymer NumSeqs
# Minimum:        1044    1056    2       0       2       1
# 2.5%-tile:      13862   23444   252     0       3       167533
# 25%-tile:       13862   23444   253     0       4       1675326
# Median:         13862   23444   253     0       4       3350651
# 75%-tile:       13862   23444   253     0       4       5025976
# 97.5%-tile:     13862   23444   253     0       6       6533768
# Maximum:        43115   43116   276     0       51      6701300
# Mean:           13866.7 23443.3 252.802 0       4.15174
# of Seqs:        6701300

# How to interpret: Most of the sequences start at position 13’862 and end at position 23,444. 
# Some sequences in the example start at position 1044 or 43115 and end at 1056 or 43116.
# These deviants from the mode positions are likely due to an insertion or deletion at the 
# terminal ends of the alignments. Sometimes you'll see sequences that start and end at the
# same position indicating a very poor alignment, which is generally due to non-specific
# amplification. 

# C) Run screen.seqs:
# To make sure that everything overlaps the same region we'll run screen.seqs to get sequences
# that start at or before position 13862 and end at or after position 23444 (which is based 
# on the summary results). We'll also set the maximum homopolymer length to 8 since there's 
# nothing in the database with a stretch of 9 or more of the same base in a row.
# From the example above:

screen.seqs(fasta=trimmed_seqs.align, start=13862, end=23444, maxhomop=8, processors=24)

Output: trimmed_seqs.good.align

# D) Make sure that our sequences only overlap the specific region:
# We filter the sequences to remove the overhangs at both ends. Since we've done paired-end
# sequencing, this shouldn't be much of an issue. In addition, there are many columns in 
# the alignment that only contain gap characters (i.e. "-"). These can be pulled out without 
# losing any information. We'll do all this with:

filter.seqs(fasta=trimmed_seqs.good.align, vertical=T, trump=.)

Output: trimmed_seqs.good.filter.fasta

# Edit the file to convert the alignment to fasta file (Get rid of the “-”) using perl 
# out of MOTHUR

perl -pe 's/-//g' trimmed_seqs.good.filter.fasta > trimmed_seqs_MOTHUR.fna

## 8.- OTU picking or clustering of sequences into OTUs using UPARSE from USEARCH at 97% 
#  sequence identity 

# A) Deduplicate the sequences (it is like clustering at 100% sequence identity, but it 
# keeps the abundance information in the header).

usearch -derep_fulllength trimmed_seqs_MOTHUR.fna -fastaout trimmed_seqs_MOTHUR_unique.fna \
 -sizeout -minseqlength 64 -threads 20
 
# Note:  When you have to many samples, usearch won’t be able to run this step, you can use 
# VSEARCH instead and then use Usearch or vsearch for the following steps:

vsearch --derep_full trimmed_seqs_MOTHUR.fna --output trimmed_seqs_MOTHUR_unique_VSEARCH.fna \
 --log=log --sizeout --minseqlength 64
 
# B) Make the reference (chimeras are filtered in this step):
# Options: -minsize 2 will removed singletons, if its 3 it will remove doubletons and so on
usearch -cluster_otus trimmed_seqs_MOTHUR_unique.fna -minsize 2 -otus otus_mc2.fa -relabel Otu

# C) Edit the header of the original file (non-deuniqued) to a format suitable for usearch:
perl -pe 'if($_=~/>.+(_\d+)/) {$_=~s/(_\d+)//g}' trimmed_seqs_MOTHUR.fna  \
> trimmed_seqs_MOTHUR_renamed.fna

# D) Make OTUs:
## Assign sequences to OTUs, with 97% cut-off (-id).
usearch -usearch_global trimmed_seqs_MOTHUR_renamed.fna -db otus_mc2.fa \
-strand plus -id 0.97 -otutabout otutab_mc2.txt 

## 9.- Assign taxonomy to the uniqued sequences
# The Greengenes database of 16S sequences is the database of reference 16S sequences used
# to assign the taxonomy. A Qiime python script is used for this with the file 97_otus.fasta
# that functions as a reference FASTA file of all sequences with known taxonomy.
parallel_assign_taxonomy_uclust.py -i otus_mc2.fa -o tax_otus_mc2 -O 20 \
 -t $PATH/qiime-deploy/qiime_software/gg_otus-13_8-release/taxonomy/97_otu_taxonomy.txt \
 -r $PATH/qiime-deploy/qiime_software/gg_otus-13_8-release/rep_set/97_otus.fasta
 
# Output in tax_otus_mc2: 
# otus_mc2_tax_assignments.log  otus_mc2_tax_assignments.txt

## 10.-	Make a biom file and add the taxonomic and metadata information:

# A) Convert a tab-delimited table to a JSON biom format. The biom format is designed to
# be a general-use format for representing biological sample by observation contingency tables 
biom convert -i otutab_mc2.txt -o otutab_mc2.biom --table-type="OTU table" --to-json

# B) Add the specific header to file with taxonomies (the # symbol has to be in the file)
# and word are separated by tabs:
nano head
#OTUID  taxonomy        confidence

cat head  otus_mc2_tax_assignments.txt > otus_mc2_tax_assignments_C.txt

# C) Add the taxonomy to the biom file:
biom add-metadata --sc-separated taxonomy --observation-header OTUID,taxonomy \
--observation-metadata-fp $PATH/otus_mc2_tax_assignments_C.txt -i otutab_mc2.biom \
-o otutab_mc2_tax.biom

# D) Add metadata:
biom add-metadata -i otutab_mc2_tax.biom -o otutab_mc2_AllMeta.biom \
--sample-metadata-fp Mapping_file.csv

# Note:  if you have an error it would probably be related to not having the same samples
# in the sequences and in the mapping file. Fix that and repeat

## 11.- Obtain a phylogenetic tree: 
# A) Align sequences from the reference using Mothur

mothur

align.seqs(fasta=otus_mc2.fa, reference=silva.bacteria.fasta, flip=t, processors=24)

#Output:  otus_mc2.align

# B) Make sure that your sequences only overlap the specific region:
filter.seqs(fasta=otus_mc2.align, vertical=T, trump=.)

# Output: otus_mc2.filter.fasta

# C) Make the phylogenetic tree using qiime (out of MOTHUR):
make_phylogeny.py -i otus_mc2.filter.fasta -o fasttree_mc2

# Note: The phylogenetic tree is added in phyloseq

# Output: fasttree_mc2.tre

## 12.- For making ordination plots you will need to normalize the .biom table. This can be
#  done using the CSS normalization in QIIME (from “Robust methods for differential abundance
#  analysis in marker genes surveys” to do PCoA plots):
normalize_table.py -i otutab_mc2_AllMeta.biom -a CSS –o CSS_normalized_otutab_mc2_AllMeta.biom

# Note: This step is not necessary for differential abundance (DA) analysis. For this the 
# un-normalized table is used and probably the DA program will do the normalization