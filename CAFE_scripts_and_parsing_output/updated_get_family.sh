#!/bin/bash
# requirements

# 1: you need the protein master file of all the sequences
# that is set here
proteins=new_genomes_isofilt_startstop_filt_blast_input.fa
# 2: Before running this script, you need to make an index of the master fasta file (with all the protein sequences)
# to do this, run, for example"
# samtools faidx makeblastdb_input.fa

# get_family.sh
# 1 = family number
# 2 = MCL dump file

# get the line number, change tabs to returns, sort, save as list
sed ''$1'q;d' $2 | tr "\t" "\n" | sort > $1.sorted_list
# now use the list to grab proteins

# search and extract
# assuming that the master fasta file is "makeblastdb_input.fa".
# if you using a different protein list, you will need to replace that name
# above, where proteins=new_fasta_file. Don't for get to index it
while read p; do samtools faidx $proteins $p ; done < $1.sorted_list > $1.family.$proteins

# run example below
# ./updated_get_family.sh 206 ../dump.new_genomes_Aug23_blast_output.mci.I30
