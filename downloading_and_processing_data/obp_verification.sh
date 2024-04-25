# Just as before, we want to get information on lengths of all the genes in a family
# Family 264 in this case, the OBPs.

cd /mnt/griffin/handor/2023_CAFE_redo/all_v_all_blast/output/node13_getfams/butterflyrate_node13_sigs
mkdir OBP_check
cd OBP_check
ln -s /mnt/griffin/handor/264.family.new_genomes_isofilt_startstop_filt_blast_input.fa .
/data/programs/exonerate-2.2.0-x86_64/bin/fastalength 264.family.new_genomes_isofilt_startstop_filt_blast_input.fa > 264.lengths
more *lengths

sed 's/G000.*.$//g' 264.lengths | awk '{gsub(" ","\t",$0); print;}' > 264.lengths.tsv
file_to_head=264.lengths.tsv
echo 'length,species' |awk '{gsub(",","\t",$0); print;}' >> $file_to_head
sed -i '1h;1d;$!H;$!d;G' $file_to_head
# download for plotting

mkdir validation264
cd validation264

ln -s /mnt/griffin/handor/2023_CAFE_redo/updated_genomes/Cydia_splendana-GCA_910591565.1.fa Cydspl.fa
ln -s /mnt/griffin/handor/2023_CAFE_redo/updated_genomes/Pluxyl_3.1_genome.fa Pluxyl.fa
ln -s /mnt/griffin/handor/2023_CAFE_redo/updated_genomes/Spofru_corn_v7_genome.fa Spofru.fa
ln -s /mnt/griffin/handor/2023_CAFE_redo/updated_genomes/Junonia_coenia_Jc_v2.scaffolds.fa Juncoe.fa

ln -s /mnt/griffin/handor/Lep_Evolution/NCBI_genomes/Ostrina_nubilalis/GCA_008921685.1_ASM892168v1_genomic.fna Ostnub.fa
ln -s /mnt/griffin/handor/Lep_Evolution/NCBI_genomes/Trichoplusia_ni/GCF_003590095.1_tn1_genomic.fna Tricni.fa
ln -s /mnt/griffin/handor/Lep_Evolution/NCBI_genomes/Arctia_plantaginis/GCA_902825455.1_iArcPla.TrioY.curated.20190705_genomic.fna Arcpla.fa
ln -s /mnt/griffin/handor/Lep_Evolution/NCBI_genomes/Adoxophyes_honmai/GCA_005406045.1_Ah_1.0_genomic.fna Adohon.fa
ln -s /mnt/griffin/handor/Lep_Evolution/NCBI_genomes/Hyles_vespertilio/GCA_009982885.1_ASM998288v1_genomic.fna Hylves.fa
ln -s /mnt/griffin/handor/Lep_Evolution/NCBI_genomes/Hyposmocoma_kahamanoa/GCF_003589595.1_ASM358959v1_genomic.fna Hypkah.fa
ln -s /mnt/griffin/handor/Lep_Evolution/LepBase/Lep_genomes/Amyelois_transitella_v1_-_scaffolds.fa Amytra.fa
ln -s /mnt/griffin/handor/Lep_Evolution/LepBase/Lep_genomes/Bicyclus_anynana_v1.2_-_scaffolds.fa Bicany.fa
ln -s /mnt/griffin/handor/Lep_Evolution/LepBase/Lep_genomes/Bombyx_mori_ASM15162v1_-_scaffolds.fa Bommor.fa
ln -s /mnt/griffin/handor/Lep_Evolution/LepBase/Lep_genomes/Calycopis_cecrops_v1.1_-_scaffolds.fa Calcec.fa
ln -s /mnt/griffin/handor/Lep_Evolution/LepBase/Lep_genomes/Danaus_plexippus_v3_-_scaffolds.fa Danple.fa
ln -s /mnt/griffin/handor/Lep_Evolution/LepBase/Lep_genomes/Heliconius_erato_demophoon_v1_-_scaffolds.fa Heldem.fa
ln -s /mnt/griffin/handor/Lep_Evolution/LepBase/Lep_genomes/Heliconius_erato_lativitta_v1_-_scaffolds.fa Hellat.fa
ln -s /mnt/griffin/handor/Lep_Evolution/LepBase/Lep_genomes/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa Helmel.fa
ln -s /mnt/griffin/handor/Lep_Evolution/LepBase/Lep_genomes/Lerema_accius_v1.1_-_scaffolds.fa Leracc.fa
ln -s /mnt/griffin/handor/Lep_Evolution/LepBase/Lep_genomes/Manduca_sexta_Msex_1.0_-_scaffolds.fa Mansex.fa
ln -s /mnt/griffin/handor/Lep_Evolution/LepBase/Lep_genomes/Melitaea_cinxia_-_scaffolds.fa Melcin.fa
ln -s /mnt/griffin/handor/Lep_Evolution/LepBase/Lep_genomes/Operophtera_brumata_v1_-_scaffolds.fa Opebru.fa
ln -s /mnt/griffin/handor/Lep_Evolution/LepBase/Lep_genomes/Papilio_glaucus_v1.1_-_scaffolds.fa Papgla.fa
ln -s /mnt/griffin/handor/Lep_Evolution/LepBase/Lep_genomes/Papilio_machaon_Pap_ma_1.0_-_scaffolds.fa Papmac.fa
ln -s /mnt/griffin/handor/Lep_Evolution/LepBase/Lep_genomes/Papilio_polytes_Ppol_1.0_-_scaffolds.fa Pappol.fa
ln -s /mnt/griffin/handor/Lep_Evolution/LepBase/Lep_genomes/Papilio_xuthus_Pxut_1.0_-_scaffolds.fa Papxut.fa
ln -s /mnt/griffin/handor/Lep_Evolution/LepBase/Lep_genomes/Phoebis_sennae_v1.1_-_scaffolds.fa Phosen.fa
ln -s /mnt/griffin/handor/Lep_Evolution/LepBase/Lep_genomes/Pieris_napi_v1.1.scaffolds.fa Pienap.fa
ln -s /mnt/griffin/handor/Lep_Evolution/LepBase/Lep_genomes/Plodia_interpunctella_v1_-_scaffolds.fa Ploint.fa
ln -s /mnt/griffin/handor/Lep_Evolution/LepBase/Vc_v1/Vcar.v1.scaf.fa Vancar.fa

# make databases
# /data/programs/diamond_v0.9.10/diamond makedb --in GenesetA.BGIBMGAnamed.90.fa -d GenesetA.BGIBMGAnamed.90
parallel '/data/programs/ncbi-blast-2.2.31+/bin/makeblastdb -in {} -dbtype nucl' ::: *.fa
# fasta of all OBPs identified
cat ../264.family.new_genomes_isofilt_startstop_filt_blast_input.fa > 264.fna
# blast away
ls *.fa | parallel '/data/programs/ncbi-blast-2.2.31+/bin/tblastn -query 264.fna -db {} -evalue .00001 -out 264_v_{.}.tsv -outfmt 6'

# check blast output for inferred number of loci
# 3 = identity, 4 is alignment length, 11 is evalue, 12 is bit, 5 mistmatches, 6 gap_openings

blast_filter(){
	awk -F"\t" '{if ($3 > 60 && $4 > 65 && $11 < 1e-10 && $12 > 80) print $0}' "$1" |  awk '{if($1!=temp){print $0;temp=$1}}'
}
export -f blast_filter
parallel 'blast_filter {} > {.}.filtered.besthits' ::: *.tsv


count_rows(){
	wc -l "$1" | sed 's/264_v_//g' | awk ' {gsub(".filtered.besthits","",$0); print $2,$1 ;}'
}
export count_rows
for file in *.besthits; do count_rows $file; done

Adohon 18
Amytra 28
Arcpla 21
Bicany 20
Bommor 17
Calcec 6
Cydspl 14
Danple 9
Heldem 31
Hellat 30
Helmel 30
Hylves 27
Hypkah 10
Juncoe 26
Leracc 7
Mansex 21
Melcin 25
Opebru 27
Ostnub 19
Papgla 17
Papmac 21
Pappol 25
Papxut 21
Phosen 17
Pienap 19
Ploint 27
Pluxyl 22
Spofru 23
Tricni 15
Vancar 17

# Some FASTA sequences for quick NCBI BLASTs

>AdohonG00000011758.1
MYDISKDFLANNSVLPGTLPNVIPFIGKALSTIQKTLVQAQFLSKGLICVKNNPLTLEDI
NTFKMLKMPEGDHAKCFAACLFKNIGILDDMGKLTSSGARQSAKQVFANDESSLSKVENI
VQECSKVNDEEVKDGDKGCERAALAFACLTQVGPKYGLDLQF
>AdohonG00000011759.1
MTDDQKAIIRQHFEQLGMECMKDAPITSDDVANLRAKKVPSGPNAPCFLACILKKSGIMD
QSGMLQKETILEKAKQVFDDEEELKSIESYLHSCSHINNEAVSDGEKGCERAIMSYKCML
DNAAQFGFDV
>OpebruG00000014903.1
MECNKQHPISPEEMMMMKDHKMPDSENAKCLMACVMRKATFIDDKGNFSVENAIAWAGKE
FQDEPKRLESSKSIYDICKKVNDEPVSDGEKGCDRAFLLSKCLIENAPKVTLLSRYPYNT
>OpebruG00000019331.1
MSQLISVCFVLFVCCVVKYEAAITEAQKSKIHTKLLTSGLDCIKDHPLSFEHIRAFRERK
VPEDEEAKCFVKCIYKNLGIMEDNGKLSEKKARESARKIFKEGDDMLASVDQIIDECIKV
NDEETSDGEKGCDRAKLAFQCFVKHAPKLGLDVDF
>PapxutG00000017047.1
MSFVKLVICFSFIGGALARTEAEVKEFFLKQSVECTKDHPVTAEEMTMLNKHELPDSKNA
RCLLACVYRRTTWMDEKGMFMKENAYKLAQEKHPEDKAKLEKSKELFELCSKVNEETFSD
GEEGCERAAKLTKCLTENAPKMGFELD
>PapxutG00000017053.1
MSKLYCVFLFLGLAVSLRHVRALSQEDIAAIKTGLRPLIAECGKEFGVDEADIKKAKESG
KIESLDPCLFACIGKKMGMINDKGEFDVEKSSETVKKFVTDKDEQKQILEIIEKCSSVND
EAVSDDKGCDRAVLLHKCMEPYKDQFDFSK
