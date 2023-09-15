# get your genome
ln -s ../../NCBI_genomes/name/file.

# look for hidden characters and blank lines in the contig names
grep '>' -B2 -A2  Amyelois_transitella_v1_-_scaffolds.fa | sed -n 'l' | tail

# notice the '\' at the end of the fasta header, it is a hidden new line character

# remove hidden new line char from header as well as all the metadata post the scaffold name
sed 's/ .*//g' Amyelois_transitella_v1_-_scaffolds.fa > Amytra_no_ws.fasta

# define a path for the Red software, which softmasks the genome
export PATH=/data/programs/RED/redUnix64/:$PATH
# run the script with your genome and output header named (this output header will be used for names of the new files generated).
python /data/programs/scripts/redmask.py -i Amytra_no_ws.fasta -o Amytra_nows_RED > results_amytra.log

# see what you got
more results_amytra.log

[Oct 16 02:32 PM] Running Python v2.7.16
[Oct 16 02:33 PM] Loading assembly:
    Contigs: 7,301
    Length:  406,468,287 bp
    nN50:    1,586,980 bp
[Oct 16 02:33 PM] Splitting genome assembly into training set (contigs > 1000 bp)
[Oct 16 02:33 PM] Finding repeats with Red (REpeat Detector)
[Oct 16 02:43 PM] Collecting results from Red
[Oct 16 02:43 PM] Summarizing results and converting to BED format

Masked genome: /mnt/griffin/handor/braker_work/masking/Amytra_nows_RED.softmasked.fa
Repeat BED file: /mnt/griffin/handor/braker_work/masking/Amytra_nows_RED.repeats.bed
Repeat FASTA file: /mnt/griffin/handor/braker_work/masking/Amytra_nows_RED.repeats.fasta
num scaffolds: 7,301
assembly size: 406,468,287 bp
masked repeats: 128,448,884 bp (31.60%)

# make a new folder for braker work
# bring in a softmasked genome that is ready for annotation, via a soft link
ln -s /mnt/griffin/handor/braker_work/masking/Amytra_nows_RED.softmasked.fa

###########
# define your genome and protein sets.
genome=Amytra_nows_RED.softmasked.fa
proteins=/mnt/griffin/handor/braker_work/braker2_test/odb10_arthropoda_proteins.fa
###########
# now set up the other dependencies.
# install a local "key" for running the gmes Gene Mark software in your local home folder
# get the key
cp /data/programs/scripts/gm_key_64 gm_key
# put it in your home folder using this script
cp gm_key ~/.gm_key

###########
# the following paths will significantly change depending on your version of braker2 and recent updates
# but something like this will run braker
###########

# first this makes an augustus config file in your local folder that augustus needs to be able to write to
your_working_dir=$(pwd)
cp -r /data/programs/Augustus_v3.3.3/config/ $your_working_dir/augustus_config
export AUGUSTUS_CONFIG_PATH=$your_working_dir/augustus_config
chmod -R 777 $your_working_dir/augustus_config

# now set paths to the other tools
# export PATH=/mnt/griffin/chrwhe/software/BRAKER2_v2.1.5/scripts/:$PATH
export PATH=/data/programs/BRAKER2_v2.1.5/scripts/:$PATH
export AUGUSTUS_BIN_PATH=/data/programs/Augustus_v3.3.3/bin
export AUGUSTUS_SCRIPTS_PATH=/data/programs/Augustus_v3.3.3/scripts
export DIAMOND_PATH=/data/programs/diamond_v0.9.24/
export GENEMARK_PATH=/data/programs/gmes_linux_64.4.61_lic/
export BAMTOOLS_PATH=/data/programs/bamtools-2.5.1/bin/
export PROTHINT_PATH=/data/programs/ProtHint/bin/
export ALIGNMENT_TOOL_PATH=/data/programs/gth-1.7.0-Linux_x86_64-64bit/bin
export SAMTOOLS_PATH=/data/programs/samtools-1.10/
export MAKEHUB_PATH=/data/programs/MakeHub/

braker.pl --genome=$genome --prot_seq=$proteins --softmasking --cores=30

# finished, took about 10 hours

# make a fasta from the gff
gff_file=braker/braker.gtf
reference=Amytra_nows_RED.softmasked.fa
cds_outfile=Amytra_HND.CDS.fa
prot_outfile=Amytra_HND.prot.fa
/data/programs/cufflinks-2.2.1.Linux_x86_64/gffread "$gff_file" -g "$reference" -x "$cds_outfile" -y "$prot_outfile"

# have a look and validate the gtf, do they have nice start, stop and no internal stops?
more Amytra_HND.prot.fa

# get count of genes
grep -c '>' Amytra_HND.prot.fa
# for only one isoform per gene
grep '>' Amytra_HND.prot.fa | grep '.t1' | wc -l
