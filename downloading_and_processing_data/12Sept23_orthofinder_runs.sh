#### Orthofinder runs #####

# Sometimes, you might want to run your protein sets through a tool called Orthofinder
# As it can give you a LOT of useful statistics.

# Here, we do orthofinder runs for our filtered braker2 protein sets
# And our raw native protein sets.

# A nice tutorial to you get you started can be found:
# https://davidemms.github.io/orthofinder_tutorials/running-an-example-orthofinder-analysis.html

# First, let's do the filtered braker2 set, since that's the cleanest
# And thus the most straightforward.

# Move to a folder with all the protein files.
cd /mnt/griffin/handor/2023_CAFE_redo/all_v_all_blast/protein_links/just_links

# Set paths to the software
export PATH=/mnt/griffin/chrwhe/software/OrthoFinder/bin/:$PATH
orthofinder_path=/mnt/griffin/chrwhe/software/OrthoFinder

# Run on 20 cores.
$orthofinder_path/orthofinder -f /mnt/griffin/handor/2023_CAFE_redo/all_v_all_blast/protein_links/just_links -a 20

# Now, getting the raw native sets in to orthofinder is a bit difficult
# Since you can't have any periods in the file.
# Or weird spaces anywhere.

# So, we need a script that removes all terminal stop codons (if they exist)
# as well as cuts out sequences with internal stops
# And removes spaces.

nano startnostopnospace.sh

#!/bin/bash
FILES="/mnt/griffin/handor/2023_CAFE_redo/updated_genomes/native_protein_sets/orthofinder_prep/*.fa*"
for f in $FILES
do
  awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $f | \
  awk '{gsub(/\.$/,"",$0); print;}' | \
  awk -F '\t'  '($2 !~ /\./)' |\
  tr "\t" "\n"| tr -d ' '>"$f.startnostop.nospace.fa"
done

./startnostopnospace.sh

# Now move to the folder you want to run OF from.
cd /mnt/griffin/handor/2023_CAFE_redo/raw_orthofinder

# Link the lightly reformatted proteins.
ln -s /mnt/griffin/handor/2023_CAFE_redo/updated_genomes/native_protein_sets/orthofinder_prep/*.startnostop.nospace.fa .

# Now we should be able to run things.
export PATH=/mnt/griffin/chrwhe/software/OrthoFinder/bin/:$PATH
orthofinder_path=/mnt/griffin/chrwhe/software/OrthoFinder
# OrthoFinder version 2.5.1 Copyright (C) 2014 David Emms
$orthofinder_path/orthofinder -f /mnt/griffin/handor/2023_CAFE_redo/raw_orthofinder -a 20
