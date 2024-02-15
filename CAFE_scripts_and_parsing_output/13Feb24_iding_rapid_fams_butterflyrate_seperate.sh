# Now that we have output files for our multi-lambda runs
# Let's see what the rapid families are at Node 13 (butterfly crown)

# Go into our cafe_reports folder.
cd /mnt/griffin/handor/2023_CAFE_redo/cafe_runs/cafe_reports/butterflyrate

# What we're interested in parsing are the lambdamu files ending in *.cafe.
# Using a script by Gregg Thomas, one of CAFE's authors.
# I'm going to use lambdamu_r3.cafe as my input here, as that had the median lambda

# Generate reports for ALL families with
python3 /mnt/griffin/handor/Lep_Evolution/report_analysis.py -i butterflyrate_lambdamu_r3.cafe -p butterflyrate_lambdamu_r3_full

# Reports for rapid families only -- useful for what we're interested in!
python3 /mnt/griffin/handor/Lep_Evolution/report_analysis.py -i butterflyrate_lambdamu_r3.cafe -p butterflyrate_rapids_lambdamu_r3 --rapid

# So, what are the rapid families at node 13?
# Our new rapid fams
# <13>:   39[+3*],49[+2*],76[+3*],79[+3*],168[+3*],280[+2*],11[-3*],41[-3*],46[-6*],66[-6*],94[-2*],128[-4*],138[-2*],173[-2*],227[-2*],264[-2*],358[-2*],369[-2*],489[-2*],492[-2*],527[-2*],3077[-1*]
# Our rapids from the previous, one lambda run, for comp.
# <13>:   49[+2*],123[+2*],168[+3*],11[-3*],41[-3*],46[-7*],66[-6*],78[-3*],84[-2*],94[-2*],128[-4*],138[-2*],173[-2*],227[-2*],358[-2*],369[-2*],489[-2*],492[-2*],527[-2*],3077[-1*]

# The new fams we need to look up are 39, 76, 79, 280, and 264
# Let's do it!

# Here, we use a script written for exploration of the single lambda data

cd /mnt/griffin/handor/2023_CAFE_redo/all_v_all_blast/output/node13_getfams
./updated_get_family.sh 39 ../dump.new_genomes_Aug23_blast_output.mci.I30
./updated_get_family.sh 76 ../dump.new_genomes_Aug23_blast_output.mci.I30
./updated_get_family.sh 79 ../dump.new_genomes_Aug23_blast_output.mci.I30
./updated_get_family.sh 280 ../dump.new_genomes_Aug23_blast_output.mci.I30
./updated_get_family.sh 264 ../dump.new_genomes_Aug23_blast_output.mci.I30

# If we want to know the different Viterbi p values for families at node 13, what can we do?
# We'll need to parse the .cafe file.
grep '^280' butterflyrate_lambdamu_r3.cafe|cut -f 1,4 

# The first number in the 19th set of parentheses is our p-value for node 13.

# What if we want counts at each node?
# Family ID       Pluxyl<0>       Adohon<2>       Cydspl<4>       <3>     Papgla<6>       Pappol<8>       Papmac<10>      Papxut<12>    <11>    <9>     <7>     Leracc<14>      Phosen<16>      Pienap<18>      <17>    Calcec<20>      Danple<22>      Bicany<24>    Helmel<26>      Heldem<28>      Hellat<30>      <29>    <27>    Juncoe<32>      Melcin<34>      <33>    Vancar<36>   <35>     <31>    <25>    <23>    <21>    <19>    <15>    <13>    Hypkah<38>      Ostnub<40>      Amytra<42>      Ploint<44>   <43>     <41>    Bommor<46>      Hylves<48>      Mansex<50>      <49>    <47>    Opebru<52>      <51>    Arcpla<54>      Spofru<56>    Tricni<58>      <57>    <55>    <53>    <45>    <39>    <37>    <5>

# 264     9       5       7       6       5       4       5       3       4       4       4       4       5       4       4    64       4       4       4       6       4       4       3       5       4       4       4       4       4       4       4    44       4       9       8       7       6       6       6       8       5       4       5       6       7       6       9    18       8       8       8       6       6       6       6       6
# 489   14      11      15      10      1       1       1       2       1       1       1       0       1       1       1    01       1       3       1       1       1       1       5       0       1       2       1       1       1       1       1    11       1       2       4       9       5       5       3       7       2       2       2       3       2       3       3    62       3       3       3       3       3       3       4
