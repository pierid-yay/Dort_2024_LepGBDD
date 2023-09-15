# 09 Sept 2023
# Parsing cafe output for downstream analyses

# Now that the cafe runs are done, it's time to put the output in our Bayesian model
# But first we need to parse it and pair it with our host order counts.

#  A tsv with those counts is located:
/mnt/griffin/handor/2023_CAFE_redo/Sept2023_host_orders.tsv

head Sept2023_host_orders.tsv
# Family_ID       Num_Orders
# Pluxyl  1
# Adohon  6
# Cydspl  1
# Papgla  4
# Pappol  1
# Papmac  3
# Papxut  1
# Leracc  1
# Phosen  1

# Great, time to get parsing!
cd /mnt/griffin/handor/2023_CAFE_redo/cafe_runs/cafe_reports

# We can generate fulll reports on counts at each node with a command like:
python3 /mnt/griffin/handor/Lep_Evolution/report_analysis.py -i lambdamu_r1.cafe -p lambdamu_r1_report

# The report we're going to parse today is a .count file for a lambdamu run.
file=lambdamu_r1_report.count

cat $file | \
# get rid of space in Family ID, make everything column separated
# then transpose using rs (reshape) command, then take only tip data
sed 's/Family ID/Family_ID/g' | column -t | rs -T | grep -v '^<' |\
# then get rid of tip numbers, squeeze all adjacent spaces into just one,
sed 's/<.*.>//g' | tr -s " " | \
# change single space to tab
tr " " "\t"  > counts_new_genomes_lambdamu_p05_r1.tsv

# using lookup table to replace via awk
# then when file2 is being read, we create an array "a" which stores second column
# from file2, indexed based upon the first column. This is like a dictionary in python
# The we print file1 row in its entirety ($0) followed by indicating which column of file1
# We want to have looked up in the dictionary and replaced.  Here we are using
# column 1 of file1, which we indcate as a[$1].
# Also, we're adding family headers

file2=/mnt/griffin/handor/2023_CAFE_redo/Sept2023_host_orders.tsv
file1=counts_new_genomes_lambdamu_p05_r1.tsv
awk 'BEGIN {FS=OFS="\t"}  NR==FNR{a[$1]=$2;next}{print $0,a[$1]}' $file2 $file1 | \
awk 'BEGIN {FS=OFS="\t"} {if($1~/^Family_ID/){gsub(/\t/,"\tFam_",$0);print $0}else{;print $0}}' |\
sed 's/Fam_Num_Orders/Num_Orders/g' | column -t > orders_appended_counts_new_genomes_lambdamu_p05_r1.tsv
sed 's/Family_ID/Species/g' orders_appended_counts_new_genomes_lambdamu_p05_r1.tsv > R_ready_orders_appended_counts_new_genomes_lambdamu_p05_r1.tsv

# Look at the file a bit
nano counts_AGATiso_b2_lambmu_error_filt_p05_r2_FOR_R.tsv
# Looks nice!



###### For large fams included ########
cd /mnt/griffin/handor/2023_CAFE_redo/cafe_runs/cafe_reports

# Generate fulll reports on counts at each node:
python3 /mnt/griffin/handor/Lep_Evolution/report_analysis.py -i largefams_lambdamu_r1_corrected.cafe -p largefams_lambdamu_r1_report

# The report we're going to parse today is a .count file for the large run.
file=largefams_lambdamu_r1_report.count

cat $file | \
# get rid of space in Family ID, make everything column separated
# then transpose using rs (reshape) command, then take only tip data
sed 's/Family ID/Family_ID/g' | column -t | rs -T | grep -v '^<' |\
# then get rid of tip numbers, squeeze all adjacent spaces into just one,
sed 's/<.*.>//g' | tr -s " " | \
# change single space to tab
tr " " "\t"  > counts_largefams_lambdamu_p05_r1.tsv

# using lookup table to replace via awk
# then when file2 is being read, we create an array "a" which stores second column
# from file2, indexed based upon the first column. This is like a dictionary in python
# The we print file1 row in its entirety ($0) followed by indicating which column of file1
# We want to have looked up in the dictionary and replaced.  Here we are using
# column 1 of file1, which we indcate as a[$1].
# Also, we're adding family headers

file2=/mnt/griffin/handor/2023_CAFE_redo/Sept2023_host_orders.tsv
file1=counts_largefams_lambdamu_p05_r1.tsv
awk 'BEGIN {FS=OFS="\t"}  NR==FNR{a[$1]=$2;next}{print $0,a[$1]}' $file2 $file1 | \
awk 'BEGIN {FS=OFS="\t"} {if($1~/^Family_ID/){gsub(/\t/,"\tFam_",$0);print $0}else{;print $0}}' |\
sed 's/Fam_Num_Orders/Num_Orders/g' | column -t > orders_appended_counts_largefams_lambdamu_p05_r1.tsv
sed 's/Family_ID/Species/g' orders_appended_counts_largefams_lambdamu_p05_r1.tsv > R_ready_orders_appended_counts_largefams_lambdamu_p05_r1.tsv

# Look at the file a bit
nano R_ready_orders_appended_counts_largefams_lambdamu_p05_r1.tsv
