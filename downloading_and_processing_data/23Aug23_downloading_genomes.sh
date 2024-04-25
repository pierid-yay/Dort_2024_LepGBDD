#### Walkthrough for downloading, assessing, annotating, and filtering genomes/protein sets ###

### DOWNLOADING UPDATED GENOMES ###
# Per reviewer request, we downloaded updated genomes for four species
# Due to the high duplication rates in the genomes used for our original submission
# We could not find a high quality genome for Tricni, so we unfortunately couldn't update that.
# There was a much higher quality genome for Cydia splendana than for pomonella, so we switched species.

cd /mnt/griffin/handor/2023_CAFE_redo/updated_genomes

# Pluxyl
# NCBI page: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_932276165.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/932/276/165/GCF_932276165.1_ilPluXylo3.1/GCF_932276165.1_ilPluXylo3.1_genomic.fna.gz
gunzip *.gz
mv GCF_932276165.1_ilPluXylo3.1_genomic.fna Pluxyl_3.1_genome.fa

# Spofru
wget https://bipaa.genouest.org/sp/spodoptera_frugiperda_pub/download/genome/corn/v7.0/sfC.ver7.fa
mv sfC.ver7.fa Spofru_corn_v7_genome.fa

# Cydspl
wget https://ftp.ensembl.org/pub/rapid-release/species/Cydia_splendana/GCA_910591565.1/ensembl/genome/Cydia_splendana-GCA_910591565.1-softmasked.fa.gz
gunzip Cydia_splendana-GCA_910591565.1-softmasked.fa.gz

# Juncoe
# Here, we use the alternate assembly available on Lepbase
# http://ensembl.lepbase.org/Junonia_coenia_jcv2/Info/Index
wget http://download.lepbase.org/v4/sequence/Junonia_coenia_Jc_v2.scaffolds.fa.gz
gunzip *.gz

### NATIVE ANNOTATIONS & PROTEIN SETS ####
# When available, I downloaded annotations and/or protein sets available with each genome
# They are located in /mnt/griffin/handor/2023_CAFE_redo/updated_genomes/native_gffs
# OR /mnt/griffin/handor/2023_CAFE_redo/updated_genomes/native_proteins

wget http://download.lepbase.org/v4/sequence/Junonia_coenia_Jc_v2.proteins.fa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/932/276/165/GCF_932276165.1_ilPluXylo3.1/GCF_932276165.1_ilPluXylo3.1_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/932/276/165/GCF_932276165.1_ilPluXylo3.1/GCF_932276165.1_ilPluXylo3.1_protein.faa.gz
wget https://ftp.ensembl.org/pub/rapid-release/species/Cydia_splendana/GCA_910591565.1/ensembl/geneset/2021_11/Cydia_splendana-GCA_910591565.1-2021_11-genes.gff3.gz
wget https://ftp.ensembl.org/pub/rapid-release/species/Cydia_splendana/GCA_910591565.1/ensembl/geneset/2021_11/Cydia_splendana-GCA_910591565.1-2021_11-pep.fa.gz
wget https://bipaa.genouest.org/sp/spodoptera_frugiperda_pub/download/annotation/corn/OGS7.0_20190530/OGS7.0_20190530.gff
####

### CHECKING THE QUALITY OF GENOMES ###
# As a quick check, we can run AsmQC
/data/programs/scripts/AsmQC Cydia_splendana-GCA_910591565.1-softmasked.fa > asmqc/cydspl_asmqc_results.txt
/data/programs/scripts/AsmQC Junonia_coenia_Jc_v2.scaffolds.fa > asmqc/juncoe_asmqc_results.txt
/data/programs/scripts/AsmQC Spofru_corn_v7_genome.fa > asmqc/spofru_asmqc_results.txt
/data/programs/scripts/AsmQC Pluxyl_3.1_genome.fa > asmqc/pluxyl_asmqc_results.txt

# But we also want to get the BUSCO statistics for these.
conda activate /mnt/griffin/handor/software/conda/miniconda3/envs/busco_env

# Run BUSCO on these genomes in genome mode.
for file in *.fa*; do
busco -i $file -l lepidoptera_odb10 -m genome -o $file.genome.busco -c 20
done

conda deactivate

### MASKING THE NEW GENOMES WITH RED ###
mkdir masking
cd masking

# Install missing RED dependencies
pip install natsort

ln -s /mnt/griffin/handor/2023_CAFE_redo/updated_genomes/GCF_932276165.1_ilPluXylo3.1_genomic.nomask.fna .
export PATH=/mnt/griffin/chrwhe/software/redmask:/home/chrwhe/.local/bin:/data/programs/redUnix64:$PATH
ingenome=GCF_932276165.1_ilPluXylo3.1_genomic.nomask.fna
outgenome=Pluxyl_MASKED_3.1_NEW_genome.fa
python2.7 /mnt/griffin/chrwhe/software/redmask/redmask.py -i $ingenome -o $outgenome > Pluxyl_NEW.log

ln -s /mnt/griffin/handor/2023_CAFE_redo/updated_genomes/Spofru_corn_v7_genome.fa .
export PATH=/mnt/griffin/chrwhe/software/redmask:/home/chrwhe/.local/bin:/data/programs/redUnix64:$PATH
ingenome=Spofru_corn_v7_genome.fa
outgenome=Spofru_MASKED_corn_v7_genome
python2.7 /mnt/griffin/chrwhe/software/redmask/redmask.py -i $ingenome -o $outgenome > Spofru.log

ln -s /mnt/griffin/handor/2023_CAFE_redo/updated_genomes/Junonia_coenia_Jc_v2.scaffolds.fa  .
export PATH=/mnt/griffin/chrwhe/software/redmask:/home/chrwhe/.local/bin:/data/programs/redUnix64:$PATH
ingenome=Junonia_coenia_Jc_v2.scaffolds.fa
outgenome=Juncoe_MASKED_v2_genome
python2.7 /mnt/griffin/chrwhe/software/redmask/redmask.py -i $ingenome -o $outgenome > Juncoe.log

ln -s /mnt/griffin/handor/2023_CAFE_redo/updated_genomes/Cydia_splendana-GCA_910591565.1.fa .
export PATH=/mnt/griffin/chrwhe/software/redmask:/home/chrwhe/.local/bin:/data/programs/redUnix64:$PATH
ingenome=Cydia_splendana-GCA_910591565.1.fa
outgenome=Cydspl_MASKED_new_genome
python2.7 /mnt/griffin/chrwhe/software/redmask/redmask.py -i $ingenome -o $outgenome > Cydspl_new.log

### BRAKER ANNOTATIONS ###
cd /mnt/griffin/handor/2023_CAFE_redo/updated_genomes/braker_annotation

wget https://v100.orthodb.org/download/odb10_arthropoda_fasta.tar.gz
tar xvf odb10_arthropoda_fasta.tar.gz
cat arthropoda/Rawdata/* > odb10_arthropoda_proteins.fa

# Bring in a masked genome that is ready for annotation, via a soft link
ln -s /mnt/griffin/handor/2023_CAFE_redo/updated_genomes/masking/Cydpom_MASKED_v2_genome.fa.softmasked.fa .

###########
# Now set up the other dependencies.
# Install a local "key" for running the gmes Gene Mark software in your local home folder
# This key needs to be updated every few months, so if it isn't working, check to see if you need a refresh.

cp /data/programs/scripts/gm_key_64 gm_key
cp gm_key ~/.gm_key

# Make an augustus config file in your local folder that augustus can write to
your_working_dir=$(pwd)
cp -r /data/programs/augustus-3.3.3/config/ $your_working_dir/augustus_config
export AUGUSTUS_CONFIG_PATH=$your_working_dir/augustus_config

# Set paths to the other tools
export PATH=/data/programs/BRAKER2_v2.1.5/scripts/:$PATH
export AUGUSTUS_BIN_PATH=/data/programs/augustus-3.3.3/bin
export AUGUSTUS_SCRIPTS_PATH=/data/programs/augustus-3.3.3/scripts
export DIAMOND_PATH=/data/programs/diamond_v0.9.10/
export GENEMARK_PATH=/mnt/griffin/chrwhe/software/braker_dependencies/gmes_linux_64.4.61_lic
export BAMTOOLS_PATH=/data/programs/bamtools-2.5.1/bin/
export PROTHINT_PATH=/data/programs/ProtHint-2.5.0/bin/
export ALIGNMENT_TOOL_PATH=/data/programs/gth-1.7.0-Linux_x86_64-64bit/bin
export SAMTOOLS_PATH=/data/programs/samtools-1.9/

###########
# Define your genome and protein sets.
genome=Cydpom_MASKED_v2_genome.fa.softmasked.fa
proteins=odb10_arthropoda_proteins.fa
# Run
braker.pl --genome=$genome --prot_seq=$proteins --softmasking --cores=30


###### FILTERING THE ANNOTATIONS #######
# So we want to filter these nice new braker annotations for correct start/stops/no internal stops
# And so that there's only one isoform per gene.
# We'll use a bash script for that.

nano annotation_filtering_v2.sh

#/bin/bash
species=$1
genome=$2
braker_gtf=$3

# example
# ./batch_TE_filtering.sh Pienap Pienap_nows_RED.softmasked.fa braker.gtf

# remove unnecesary features for EBI submission
/mnt/griffin/racste/Software/gffread-0.12.6.Linux_x86_64/gffread $braker_gtf -o "$species".temp.gff --t-adopt

# some overlapping gene parent features: use agat merge function
agat_convert_sp_gxf2gxf.pl -g "$species".temp.gff -o "$species".temp.merged.gff -ml

# Give all features consistent ID
agat_sp_manage_IDs.pl --gff "$species".temp.merged.gff --prefix "$species" --ensembl --tair --type_dependent --collective -o $species.temp.merged.rename.gff

# Filter so that it's longest isoform only
agat_sp_keep_longest_isoform.pl -gff "$species".temp.merged.rename.gff -o "$species".temp.merged.rename.longest_iso.gff

# We'll need to tweak the naming just a bit more still,
# because spaces will be introduced otherwise in the protein file
cat "$species".temp.merged.rename.longest_iso.gff | awk -F '\\;geneID=' '{print $1}' > "$species".temp.merged.rename.longest_iso.edited.gff

# Call proteins using the filtered gff from before and the hardmasked fasta we just made.
reference="$genome"
filtered_gff_file="$species".temp.merged.rename.longest_iso.edited.gff
root="$species".Jfilt.isofilt
# -yJ flag only prints converted amino acid sequence of any mRNAs with no in-exon stop codons, and with existing start and stop codon.
/data/programs/cufflinks-2.2.1.Linux_x86_64/gffread "$filtered_gff_file" -g "$reference" -J -y "$root".prot.fa

# close out of nano
chmod u+x annotation_filtering_v2.sh

# activate conda env for agat
conda activate agat_env

# run it
./annotation_filtering_v2.sh Spofru spofru/Spofru_MASKED_corn_v7_genome.softmasked.ed.fa spofru/braker/braker.gtf
./annotation_filtering_v2.sh Juncoe juncoe/Juncoe_MASKED_v2_genome.softmasked.nowhite.fa juncoe/braker/braker.gtf
./annotation_filtering_v2.sh Cydspl cydspl/Cydspl_MASKED_new_genome.softmasked.ed.fa cydspl/braker/braker.gtf
./annotation_filtering_v2.sh Pluxyl pluxyl/Pluxyl_MASKED_3.1_NEW_genome.fa.softmasked.ed.fa pluxyl/braker/braker.gtf

conda deactivate

##### Get gene counts ########
cd braker_annotation/filtered_protein_sets
grep -c '>' *.fa
# Spofru.Jfilt.isofilt.prot.fa:20466
# Juncoe.Jfilt.isofilt.prot.fa:20664
# Pluxyl.Jfilt.isofilt.prot.fa:17556
# Cydspl.Jfilt.isofilt.prot.fa:22275

##### BUSCO the protein sets ######
conda activate /mnt/griffin/handor/software/conda/miniconda3/envs/busco_env

# Run BUSCO on these genomes in protein mode.
for file in *.fa*; do
busco -i $file -l lepidoptera_odb10 -m protein -o $file.prot.busco -c 20
done

# Spofru: C:96.0%[S:93.8%,D:2.2%],F:1.8%,M:2.2%,n:5286
# Juncoe: C:96.7%[S:95.4%,D:1.3%],F:1.4%,M:1.9%,n:5286
# Pluxyl: C:93.1%[S:91.9%,D:1.2%],F:1.4%,M:5.5%,n:5286 # high missing???
# Cydspl: C:96.6%[S:95.3%,D:1.3%],F:0.9%,M:2.5%,n:5286

##### Extras: Making some non-start-stop filtered protein sets for comparisons #####
cd raw_prot_sets

nano no_start_stop_annotation_filtering.sh

#/bin/bash
species=$1
genome=$2
braker_gtf=$3

# example
# ./batch_TE_filtering.sh Pienap Pienap_nows_RED.softmasked.fa braker.gtf

# remove unnecesary features for EBI submission
/mnt/griffin/racste/Software/gffread-0.12.6.Linux_x86_64/gffread $braker_gtf -o "$species".temp.gff --t-adopt

# some overlapping gene parent features: use agat merge function
agat_convert_sp_gxf2gxf.pl -g "$species".temp.gff -o "$species".temp.merged.gff -ml

# Give all features consistent ID
agat_sp_manage_IDs.pl --gff "$species".temp.merged.gff --prefix "$species" --ensembl --tair --type_dependent --collective -o $species.temp.merged.rename.gff

# Filter so that it's longest isoform only
agat_sp_keep_longest_isoform.pl -gff "$species".temp.merged.rename.gff -o "$species".temp.merged.rename.longest_iso.gff

# We'll need to tweak the naming just a bit more still,
# because spaces will be introduced otherwise in the protein file
cat "$species".temp.merged.rename.longest_iso.gff | awk -F '\\;geneID=' '{print $1}' > "$species".final.longest_iso.edited.gff

# Call proteins using the filtered gff from before and the hardmasked fasta we just made.
reference="$genome"
filtered_gff_file="$species".final.longest_iso.edited.gff
root="$species".noJfilt.isofilt

# Note there is no -J flag here.
/data/programs/cufflinks-2.2.1.Linux_x86_64/gffread "$filtered_gff_file" -g "$reference" -y "$root".prot.fa

# close nano and chmod file.

./no_start_stop_annotation_filtering.sh Spofru spofru/Spofru_MASKED_corn_v7_genome.softmasked.ed.fa spofru/braker/braker.gtf
./no_start_stop_annotation_filtering.sh Juncoe juncoe/Juncoe_MASKED_v2_genome.softmasked.nowhite.fa juncoe/braker/braker.gtf
./no_start_stop_annotation_filtering.sh Cydspl cydspl/Cydspl_MASKED_new_genome.softmasked.ed.fa cydspl/braker/braker.gtf
./no_start_stop_annotation_filtering.sh Pluxyl pluxyl/Pluxyl_MASKED_3.1_NEW_genome.fa.softmasked.ed.fa pluxyl/braker/braker.gtf


##### And filtering them with the start/stop script, rather than the J filter with gff read #####

# First, convert everything to a two-line fasta.

python2 /data/programs/scripts/multiline_to_two_line.py Spofru.noJfilt.isofilt.prot.fa > Spofru.2line.fa
python2 /data/programs/scripts/multiline_to_two_line.py Cydspl.noJfilt.isofilt.prot.fa > Cydspl.2line.fa
python2 /data/programs/scripts/multiline_to_two_line.py Pluxyl.noJfilt.isofilt.prot.fa > Pluxyl.2line.fa
python2 /data/programs/scripts/multiline_to_two_line.py Juncoe.noJfilt.isofilt.prot.fa > Juncoe.2line.fa

nano start_stop_awkscript.sh

#!/bin/bash
FILES="/mnt/griffin/handor/2023_CAFE_redo/updated_genomes/braker_annotation/raw_prot_sets/*2line.fa"
for f in $FILES
do
  awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $f | awk -F '\t'  '($2 ~ /^M/)' | tr "\t" "\n" | \
  awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' |  awk -F '\t'  '($2 ~ /\.$/)' |  tr "\t" "\n" | \
  awk '{gsub(/\.$/,"",$0); print;}' | \
  awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' |  awk -F '\t'  '($2 !~ /\./)' |\
  tr "\t" "\n" > "${f%2line.fa}.startstop.fa"
done

# close nano and chmod the script

# run the script
./start_stop_awkscript.sh

grep -c '>' *..startstop*fa

# Cydspl..startstop.fa:22275
# Juncoe..startstop.fa:20664
# Pluxyl..startstop.fa:17556
# Spofru..startstop.fa:20466

# All counts identical to the -J filtered sets
# So the -J flag does do what we hoped it did.
/data/programs/cufflinks-2.2.1.Linux_x86_64/gffread  GCF_932276165.1_ilPluXylo3.1_genomic.gff -g Pluxyl_MASKED_3.1_NEW_genome.fa.softmasked.ed.fa -y Pluxyl.RAWnative.prot.fa
/data/programs/cufflinks-2.2.1.Linux_x86_64/gffread  OGS7.0_20190530.gff -g Spofru_corn_v7_genome.fa -J -y Spofru.startstop.noiso.prot.fa
