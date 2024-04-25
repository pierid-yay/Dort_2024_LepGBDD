# reviewer asked for some leve of validation. 
# sure, below I use their comments to direct work

# what would validation be at this point?
# This additional work is fundamental to rule out gene annotation errors as the source of what the authors
# described as “surprising correlations” between gene family expansions/contractions and phenotypic variation across lepidopterans

# Manual validation of some gene duplicates and gene losses should be performed.
# General validations include checking the length of gene duplicates that appear only in one species;
# they can be due to gene fragmentation. Gene family expansions associated in diet breath can be easily tested to rule out possible gene fragmentations.

###############
# 315	serine-type endopeptidase activity
# 315serine-type.endopeptidase.activity
315.family.new_genomes_isofilt_startstop_filt_blast_input.fa
cd validation
mkdir 315
cd 315
ln -s ../315.family.new_genomes_isofilt_startstop_filt_blast_input.fa .
../exonerate-2.2.0-x86_64/bin/fastalength 315.family.new_genomes_isofilt_startstop_filt_blast_input.fa > 315.lengths

# could plot lengths as funciton of species IDs. Most species do have similar lengths, but there are some shorter versions
# when look at the fasta file here in AliView, the shorter versions are still divergent in homologous regions, they are not simply fragments being doubly counted
# while changes in length may suggest changes in function, but nevertheless, they appear to be valid copies
# thus, the overall pattern holds of more copies, of normal full legnths in all moths compared to butterflies

# cleaning output
sed 's/G000.*.$//g' 315.lengths | awk '{gsub(" ","\t",$0); print;}' > 315.lengths.tsv
file_to_head=315.lengths.tsv
echo 'length,species' |awk '{gsub(",","\t",$0); print;}' >> $file_to_head
sed -i '1h;1d;$!H;$!d;G' $file_to_head
head 315.*
# download for plotting

###############
# 140	trypsin-like serine protease
# 140trypsin-like.serine.protease
cd validation
mkdir 140
cd 140
ln -s ../140.family.new_genomes_isofilt_startstop_filt_blast_input.fa .
../exonerate-2.2.0-x86_64/bin/fastalength 140.family.new_genomes_isofilt_startstop_filt_blast_input.fa > 140.lengths
sed 's/G000.*.$//g' 140.lengths | awk '{gsub(" ","\t",$0); print;}' > 140.lengths.tsv
file_to_head=140.lengths.tsv
echo 'length,species' |awk '{gsub(",","\t",$0); print;}' >> $file_to_head
sed -i '1h;1d;$!H;$!d;G' $file_to_head
head 140.*
# download for plotting


#################
# Losses can be checked by Blasting orthologs against genomes where they are not annotated.
# This for example should be done for the 7tm gene family in Papilionoidea.

mkdir 7tm
cd 7tm
ln -s ../489.family.new_genomes_isofilt_startstop_filt_blast_input.fa .
../exonerate-2.2.0-x86_64/bin/fastalength 489.family.new_genomes_isofilt_startstop_filt_blast_input.fa > 489.lengths
more *lengths
# could plot lengths as funciton of species IDs. Most species do have similar lengths, but there are some shorter versions
# however, there shorter versions are still having amino acid sequences that are divergent from their respective regions of shared similarity
# thus, changes in length may suggest changes in function, but nevertheless, the overall pattern holds of more copies, of normal full legnths in all moths compared
# to butterflies

head 489.lengths
410 AdohonG00000003826.1
407 AdohonG00000003827.1
415 AdohonG00000017629.1

sed 's/G000.*.$//g' 489.lengths | awk '{gsub(" ","\t",$0); print;}' > 489.lengths.tsv
file_to_head=489.lengths.tsv
echo 'length,species' |awk '{gsub(",","\t",$0); print;}' >> $file_to_head
sed -i '1h;1d;$!H;$!d;G' $file_to_head
# download for plotting


#################
# blast assessment
# a large scale blast assessment of how many loci I can find in each species genome

# genomes

Genome_fasta	6_Letter_Code, commented out if linked below
# GCA_005406045.1_Ah_1.0_genomic.fna	Adohon
# Amyelois_transitella_v1_-_scaffolds.fa	Amytra
# GCA_902825455.1_iArcPla.TrioY.curated.20190705_genomic.fna	Arcpla
# Bicyclus_anynana_v1.2_-_scaffolds.fa	Bicany
# Bombyx_mori_ASM15162v1_-_scaffolds.fa	Bommor
# Calycopis_cecrops_v1.1_-_scaffolds.fa	Calcec
# Cydia_splendana-GCA_910591565.1-softmasked.fa	Cydspl
# Danaus_plexippus_v3_-_scaffolds.fa	Danple
# Heliconius_erato_demophoon_v1_-_scaffolds.fa	Heldem
# Heliconius_erato_lativitta_v1_-_scaffolds.fa	Hellat
# Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa	Helmel
# GCA_009982885.1_ASM998288v1_genomic.fna	Hylves
# GCF_003589595.1_ASM358959v1_genomic.fna	Hypkah
# Junonia_coenia_Jc_v2.scaffolds.fa.gz	Juncoe
# Lerema_accius_v1.1_-_scaffolds.fa	Leracc
# Manduca_sexta_Msex_1.0_-_scaffolds.fa	Mansex
# Melitaea_cinxia_-_scaffolds.fa	Melcin
# Operophtera_brumata_v1_-_scaffolds.fa	Opebru
# GCA_008921685.1_ASM892168v1_genomic.fna	Ostnub
# Papilio_glaucus_v1.1_-_scaffolds.fa	Papgla
# Papilio_machaon_Pap_ma_1.0_-_scaffolds.fa	Papmac
# Papilio_polytes_Ppol_1.0_-_scaffolds.fa	Pappol
# Papilio_xuthus_Pxut_1.0_-_scaffolds.fa	Papxut
# Phoebis_sennae_v1.1_-_scaffolds.fa	Phosen
# Pieris_napi_v1.1.scaffolds.fa	Pienap
# Plodia_interpunctella_v1_-_scaffolds.fa	Ploint
# GCF_932276165.1_ilPluXylo3.1_genomic.fna	Pluxyl
# sfC.ver7.fa	Spofru
# GCF_003590095.1_tn1_genomic.fna	Tricni
# Vcar.v1.scaf.fa	Vancar
#
cd validation
ln -s ../Lep_gene_birth_death/genomes/Cydia_splendana-GCA_910591565.1.fa Cydspl.fa
ln -s ../Lep_gene_birth_death/genomes/GCF_932276165.1_ilPluXylo3.1_genomic.fna Pluxyl.fa
ln -s ../Lep_gene_birth_death/genomes/sfC.ver7.fa Spofru.fa
ln -s ../updated_genomes/Junonia_coenia_Jc_v2.scaffolds.fa Juncoe.fa
ln -s ../Lep_Evolution/NCBI_genomes/Ostrina_nubilalis/GCA_008921685.1_ASM892168v1_genomic.fna Ostnub.fa
ln -s ../Lep_Evolution/NCBI_genomes/Trichoplusia_ni/GCF_003590095.1_tn1_genomic.fna Tricni.fa
ln -s ../Lep_Evolution/NCBI_genomes/Arctia_plantaginis/GCA_902825455.1_iArcPla.TrioY.curated.20190705_genomic.fna Arcpla.fa
ln -s ../Lep_Evolution/NCBI_genomes/Adoxophyes_honmai/GCA_005406045.1_Ah_1.0_genomic.fna Adohon.fa
ln -s ../Lep_Evolution/NCBI_genomes/Hyles_vespertilio/GCA_009982885.1_ASM998288v1_genomic.fna Hylves.fa
ln -s ../Lep_Evolution/NCBI_genomes/Hyposmocoma_kahamanoa/GCF_003589595.1_ASM358959v1_genomic.fna Hypkah.fa
ln -s ../Lep_Evolution/LepBase/Lep_genomes/Amyelois_transitella_v1_-_scaffolds.fa Amytra.fa
ln -s ../Lep_Evolution/LepBase/Lep_genomes/Bicyclus_anynana_v1.2_-_scaffolds.fa Bicany.fa
ln -s ../Lep_Evolution/LepBase/Lep_genomes/Bombyx_mori_ASM15162v1_-_scaffolds.fa Bommor.fa
ln -s ../Lep_Evolution/LepBase/Lep_genomes/Calycopis_cecrops_v1.1_-_scaffolds.fa Calcec.fa
ln -s ../Lep_Evolution/LepBase/Lep_genomes/Danaus_plexippus_v3_-_scaffolds.fa Danple.fa
ln -s ../Lep_Evolution/LepBase/Lep_genomes/Heliconius_erato_demophoon_v1_-_scaffolds.fa Heldem.fa
ln -s ../Lep_Evolution/LepBase/Lep_genomes/Heliconius_erato_lativitta_v1_-_scaffolds.fa Hellat.fa
ln -s ../Lep_Evolution/LepBase/Lep_genomes/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa Helmel.fa
ln -s ../Lep_Evolution/LepBase/Lep_genomes/Lerema_accius_v1.1_-_scaffolds.fa Leracc.fa
ln -s ../Lep_Evolution/LepBase/Lep_genomes/Manduca_sexta_Msex_1.0_-_scaffolds.fa Mansex.fa
ln -s ../Lep_Evolution/LepBase/Lep_genomes/Melitaea_cinxia_-_scaffolds.fa Melcin.fa
ln -s ../Lep_Evolution/LepBase/Lep_genomes/Operophtera_brumata_v1_-_scaffolds.fa Opebru.fa
ln -s ../Lep_Evolution/LepBase/Lep_genomes/Papilio_glaucus_v1.1_-_scaffolds.fa Papgla.fa
ln -s ../Lep_Evolution/LepBase/Lep_genomes/Papilio_machaon_Pap_ma_1.0_-_scaffolds.fa Papmac.fa
ln -s ../Lep_Evolution/LepBase/Lep_genomes/Papilio_polytes_Ppol_1.0_-_scaffolds.fa Pappol.fa
ln -s ../Lep_Evolution/LepBase/Lep_genomes/Papilio_xuthus_Pxut_1.0_-_scaffolds.fa Papxut.fa
ln -s ../Lep_Evolution/LepBase/Lep_genomes/Phoebis_sennae_v1.1_-_scaffolds.fa Phosen.fa
ln -s ../Lep_Evolution/LepBase/Lep_genomes/Pieris_napi_v1.1.scaffolds.fa Pienap.fa
ln -s ../Lep_Evolution/LepBase/Lep_genomes/Plodia_interpunctella_v1_-_scaffolds.fa Ploint.fa
ln -s ../Lep_Evolution/LepBase/Vc_v1/Vcar.v1.scaf.fa Vancar.fa

# make databases
# /data/programs/diamond_v0.9.10/diamond makedb --in GenesetA.BGIBMGAnamed.90.fa -d GenesetA.BGIBMGAnamed.90
parallel '/data/programs/ncbi-blast-2.2.31+/bin/makeblastdb -in {} -dbtype nucl' ::: *.fa
# fasta of all 7tm's identified
cat ../489.family.new_genomes_isofilt_startstop_filt_blast_input.fa > 489.fna
# blast away
ls *.fa | parallel '/data/programs/ncbi-blast-2.2.31+/bin/tblastn -query 489.fna -db {} -evalue .00001 -out 489_v_{.}.tsv -outfmt 6'

# check blast output for inferred number of loci
# 3 = identity, 4 is alignment length, 11 is evalue, 12 is bit, 5 mistmatches, 6 gap_openings

blast_filter(){
	awk -F"\t" '{if ($3 > 60 && $4 > 90 && $11 < 1e-10 && $12 > 80) print $0}' "$1" |  awk '{if($1!=temp){print $0;temp=$1}}'
}
export -f blast_filter
parallel 'blast_filter {} > {.}.filtered.besthits' ::: *.tsv


count_rows(){
	wc -l "$1" | sed 's/489_v_//g' | awk ' {gsub(".filtered.besthits","",$0); print $2,$1 ;}'
}
export count_rows
for file in *.besthits; do count_rows $file; done

# species and blast inferred loci, crude but its what the Rev wanted.
Adohon 12
Amytra 8
Arcpla 1
Bicany 1
Bommor 8
Calcec 0
Cydspl 13
Danple 1
Heldem 5
Hellat 3
Helmel 4
Hylves 1
Hypkah 2
Juncoe 6
Leracc 0
Mansex 2
Melcin 1
Opebru 1
Ostnub 4
Papgla 1
Papmac 5
Pappol 1
Papxut 1
Phosen 1
Pienap 1
Ploint 5
Pluxyl 11
Spofru 6
Tricni 2
Vancar 4

#
