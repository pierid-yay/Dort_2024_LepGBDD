##### 31 August - Blasting and clustering proteins ######

# Now that we have filtered protein sets for our new genome
# It's time to get everything ready for CAFE.

##### Make a master protein file and get the all-v-all blast running #######
# The protein sets from all other species are located in
# /mnt/griffin/handor/braker_work/RepeatMasker_removals/protein_fastas/isoform_filt_start_stop_prots/
# Note: Though the path says RepeatMasker, this folder was the control subset where RM was not run.

# I think the best thing to do would be just to softlink all the old protein sets to a new folder
# Then delete the links for the species we've gotten new genomes for
cd /mnt/griffin/handor/2023_CAFE_redo/all_v_all_blast/protein_links
ln -s /mnt/griffin/handor/braker_work/RepeatMasker_removals/protein_fastas/isoform_filt_start_stop_prots/*.fa .
rm Spofru.filt_gff_original_genome_2line.fa.startstop.fa
rm Cydpom.filt_gff_original_genome_2line.fa.startstop.fa
rm Pluxyl.filt_gff_original_genome_2line.fa.startstop.fa
rm Juncoe.filt_gff_original_genome_2line.fa.startstop.fa

# And then softlink over the four new protein sets.
ln -s /mnt/griffin/handor/2023_CAFE_redo/updated_genomes/braker_annotation/raw_prot_sets/*startstop.fa .

# Make a master protein set
cat *.fa > new_genomes_isofilt_startstop_filt_blast_input.fa

# How many proteins?
grep -c '>' new_genomes_isofilt_startstop_filt_blast_input.fa
# 581538

# Make a blast db
/data/programs/ncbi-blast-2.14.1+/bin/makeblastdb -in new_genomes_isofilt_startstop_filt_blast_input.fa  -dbtype prot -out blastdb

# Run the all-v-all blast in parallel.

# parallel blast
# to limit cores, use "-P 20"
# to set the block size, use "--block 130k"
# to keep the output in the same order as the input, use "-k"
# if using block size, tell parallel where to split the file, which for a fasta file would be the header, "--recstart '>'"
# for pipe, it is loading this at the location of the "-" mark, which is for the query here.
# if using multiple input files, and test the syntax of each call using "--dryrun", but if using pipe, as we are here, you can't use
# dryrun since there is nothing to see.

# Set the query file to the master protein set.
query=new_genomes_isofilt_startstop_filt_blast_input.fa

# Test with a subset of the file.
cat $query | head -100000 | parallel --block 130k -k --recstart '>' --pipe "/data/programs/ncbi-blast-2.14.1+/bin/blastp -num_threads 1 -db blastdb -query - -outfmt 7 -seg yes " > test_head_100k.txt

# Do the full run!
cat $query | parallel --block 130k -k --recstart '>' --pipe "/data/programs/ncbi-blast-2.14.1+/bin/blastp -num_threads 1 -db blastdb -query - -outfmt 7 -seg yes " > new_genomes_Aug23_blast_output.txt

###### MANY HOURS LATER -- MCL clustering ######
# Let's do some clustering to assign genes to "families"
# We'll use a program called mcl to do that

# convert blast output into ABC format for mcl.
grep -v "#" new_genomes_Aug23_blast_output.txt | cut -f 1,2,11 > new_genomes_Aug23_blast_output.abc

# in the following three commands, we create a network and a dictionary file, then we do the clustering

# create the tab file
/data/programs/mcl-14-137/bin/bin/mcxload -abc new_genomes_Aug23_blast_output.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o new_genomes_Aug23_blast_output.mci -write-tab new_genomes_Aug23_blast_output.tab

# Below, the-I (inflation) parameter determines how granular the clustering will be.
# Lower numbers mean denser clusters, but a value of 3 usually works
/data/programs/mcl-14-137/bin/bin/mcl new_genomes_Aug23_blast_output.mci -I 3 -te 25

#[mcl] jury pruning marks: <56,91,97>, out of 100
#[mcl] jury pruning synopsis: <69.9 or fair> (cf -scheme, -do log)
#[mclIO] writing <out.new_genomes_Aug23_blast_output.mci.I30>
#.......................................
#[mclIO] wrote native interchange 581005x33585 matrix with 581005 entries to stream <out.new_genomes_Aug23_blast_output.mci.I30>
#[mcl] 33585 clusters found
#[mcl] output is in out.new_genomes_Aug23_blast_output.mci.I30


# Make the dump file -- you'll need this for later scripts.
/data/programs/mcl-14-137/bin/bin/mcxdump -icl out.new_genomes_Aug23_blast_output.mci.I30 -tabr new_genomes_Aug23_blast_output.tab -o dump.new_genomes_Aug23_blast_output.mci.I30

# tabulate the number of gene copies in each gene family for each species
python /mnt/griffin/handor/Lep_Evolution/CAFE_tutorial/python_scripts/cafetutorial_mcl2rawcafe.py -i dump.new_genomes_Aug23_blast_output.mci.I30 -o not_sizefiltered_new_genomes_Aug23_cafe_input.txt -sp "Adohon Amytra Arcpla Bicany Bommor Calcec Cydspl Danple Heldem Hellat Helmel Hylves Hypkah Juncoe Leracc Mansex Melcin Opebru Ostnub Papgla Papmac Pappol Papxut Phosen Pienap Ploint Pluxyl Spofru Tricni Vancar"

# Check to see if they match up with the number of genes listed for a chosen family in the dump file.
# You can do this with, e.g.,
head -1 dump.new_genomes_Aug23_blast_output.mci.I30 | grep -o 'Adohon' | wc -l
# 103
# If we then look at out not_sizefiltered... file, we'll see a count of 103 for family 1. 

# filter out gene families where one or more species has more than 100 copies
# because these families can cause parameter estimates tgo be non-informative
python /mnt/griffin/handor/Lep_Evolution/CAFE_tutorial/python_scripts/cafetutorial_clade_and_size_filter.py -i not_sizefiltered_new_genomes_Aug23_cafe_input.txt -o filtered_new_genomes_Aug23_cafe_input.txt -s




########### Warning messages #########
Warning: [blastp] Query_99 AdohonG000000150.. : Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options
Warning: [blastp] Query_72 CydsplG000000141.. : Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options
Warning: [blastp] Query_197 CydsplG00000010.. : Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options
Warning: [blastp] Query_198 CydsplG00000010.. : Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options
Warning: [blastp] Query_38 MansexG000000056.. : Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options
Warning: [blastp] Query_55 OpebruG000000101.. : Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options
Warning: [blastp] Query_42 PapglaG000000181.. : Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options
Warning: [blastp] Query_24 PapxutG000000013.. : Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options
