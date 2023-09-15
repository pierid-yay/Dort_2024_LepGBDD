#### 1 September, building CAFE error model and doing some CAFE runs #####

# Now that we have clustered output, we want to build an error model for our CAFE runs.
# We do so with the following script:

cd cafe_runs
nano new_genomes_errormodel_cafe.sh

#!/mnt/griffin/handor/software/cafe
load -i ../all_v_all_blast/output/filtered_new_genomes_Aug23_cafe_input.txt -t 20 -l cafe_reports/filter_errormodel_building.log -filter
tree (Pluxyl:127.3,((Adohon:60,Cydspl:60):55.1,(((Papgla:23.8,(Pappol:14.5,(Papmac:10.8,Papxut:10.8):3.7):9.3):74.5,(Leracc:92.4,((Phosen:51.7,Pienap:51.7):33.2,(Calcec:78.6,(Danple:69.4,(Bicany:63.8,((Helmel:11.8,(Heldem:1,Hellat:1):10.8):43,((Juncoe:29.6,Melcin:29.6):8.2,Vancar:37.8):17):9):5.6):9.2):6.3):7.5):5.9):9.3,(Hypkah:105.2,((Ostnub:91.6,(Amytra:33,Ploint:33):58.6):9.7,(((Bommor:70.3,(Hylves:42.8,Mansex:42.8):27.5):18.4,Opebru:88.7):2.5,(Arcpla:68.5,(Spofru:60.6,Tricni:60.6):7.9):22.7):10.1):3.9):2.4):7.5):12.2)
lambda -s
report cafe_reports/building_errormodel_filter_r1

### Run with ###
python2 /mnt/griffin/handor/software/CAFE/cafe/caferror.py -i new_genomes_errormodel_cafe.sh -d cafe_reports/new_genomes_errormodel_filter_r1 -v 0 -f 1

#### TWO DAYS LATER #####

# We have our output, after 28 different runs.

# Global Error Estimation:         0.100390625
# Score with global errormodel:    302036.080145
# Lambda with global errormodel:   0.00213902937492

# So roughly 10% of our gene counts are expected to be incorrect.
# This is in line with the 0.05-0.15 range Gregg Thomas mentioned in a github issues thread
# So nothing out of the ordinary.

# Now that we have a completed error model, we can write our scripts for several runs #
# And all of these scripts should point to the correct error model file.
/mnt/griffin/handor/2023_CAFE_redo/cafe_runs/cafe_reports/new_genomes_errormodel_filter_r1/cafe_errormodel_0.100390625.txt

# First, ones where lambda = mu
cd cafe_runs/lambda_only/r1

nano newgenome_lambonly_cafe_r1.sh

#!/mnt/griffin/handor/software/cafe
tree (Pluxyl:127.3,((Adohon:60,Cydspl:60):55.1,(((Papgla:23.8,(Pappol:14.5,(Papmac:10.8,Papxut:10.8):3.7):9.3):74.5,(Leracc:92.4,((Phosen:51.7,Pienap:51.7):33.2,(Calcec:78.6,(Danple:69.4,(Bicany:63.8,((Helmel:11.8,(Heldem:1,Hellat:1):10.8):43,((Juncoe:29.6,Melcin:29.6):8.2,Vancar:37.8):17):9):5.6):9.2):6.3):7.5):5.9):9.3,(Hypkah:105.2,((Ostnub:91.6,(Amytra:33,Ploint:33):58.6):9.7,(((Bommor:70.3,(Hylves:42.8,Mansex:42.8):27.5):18.4,Opebru:88.7):2.5,(Arcpla:68.5,(Spofru:60.6,Tricni:60.6):7.9):22.7):10.1):3.9):2.4):7.5):12.2)
load -i ../../../all_v_all_blast/output/filtered_new_genomes_Aug23_cafe_input.txt -t 20 -l ../../cafe_reports/lambda_only_r1.log -p 0.05 -filter
errormodel -model /mnt/griffin/handor/2023_CAFE_redo/cafe_runs/cafe_reports/new_genomes_errormodel_filter_r1/cafe_errormodel_0.100390625.txt -all
lambda -s
report ../../cafe_reports/lambda_only_r1

# run with ./newgenome_lambonly_cafe_r1.sh
# Good to note, with the root filter applied
# The number of families is 9637

# Do another run
cd ../r2

nano newgenome_lambonly_cafe_r2.sh

#!/mnt/griffin/handor/software/cafe
tree (Pluxyl:127.3,((Adohon:60,Cydspl:60):55.1,(((Papgla:23.8,(Pappol:14.5,(Papmac:10.8,Papxut:10.8):3.7):9.3):74.5,(Leracc:92.4,((Phosen:51.7,Pienap:51.7):33.2,(Calcec:78.6,(Danple:69.4,(Bicany:63.8,((Helmel:11.8,(Heldem:1,Hellat:1):10.8):43,((Juncoe:29.6,Melcin:29.6):8.2,Vancar:37.8):17):9):5.6):9.2):6.3):7.5):5.9):9.3,(Hypkah:105.2,((Ostnub:91.6,(Amytra:33,Ploint:33):58.6):9.7,(((Bommor:70.3,(Hylves:42.8,Mansex:42.8):27.5):18.4,Opebru:88.7):2.5,(Arcpla:68.5,(Spofru:60.6,Tricni:60.6):7.9):22.7):10.1):3.9):2.4):7.5):12.2)
load -i ../../../all_v_all_blast/output/filtered_new_genomes_Aug23_cafe_input.txt -t 20 -l ../../cafe_reports/lambda_only_r2.log -p 0.05 -filter
errormodel -model /mnt/griffin/handor/2023_CAFE_redo/cafe_runs/cafe_reports/new_genomes_errormodel_filter_r1/cafe_errormodel_0.100390625.txt -all
lambda -s
report ../../cafe_reports/lambda_only_r2

# And a third.
cd ../r3

nano newgenome_lambonly_cafe_r3.sh

#!/mnt/griffin/handor/software/cafe
tree (Pluxyl:127.3,((Adohon:60,Cydspl:60):55.1,(((Papgla:23.8,(Pappol:14.5,(Papmac:10.8,Papxut:10.8):3.7):9.3):74.5,(Leracc:92.4,((Phosen:51.7,Pienap:51.7):33.2,(Calcec:78.6,(Danple:69.4,(Bicany:63.8,((Helmel:11.8,(Heldem:1,Hellat:1):10.8):43,((Juncoe:29.6,Melcin:29.6):8.2,Vancar:37.8):17):9):5.6):9.2):6.3):7.5):5.9):9.3,(Hypkah:105.2,((Ostnub:91.6,(Amytra:33,Ploint:33):58.6):9.7,(((Bommor:70.3,(Hylves:42.8,Mansex:42.8):27.5):18.4,Opebru:88.7):2.5,(Arcpla:68.5,(Spofru:60.6,Tricni:60.6):7.9):22.7):10.1):3.9):2.4):7.5):12.2)
load -i ../../../all_v_all_blast/output/filtered_new_genomes_Aug23_cafe_input.txt -t 20 -l ../../cafe_reports/lambda_only_r3.log -p 0.05 -filter
errormodel -model /mnt/griffin/handor/2023_CAFE_redo/cafe_runs/cafe_reports/new_genomes_errormodel_filter_r1/cafe_errormodel_0.100390625.txt -all
lambda -s
report ../../cafe_reports/lambda_only_r3

# Alright, great that those are done.
# Now it's time for a more complex model
# Where we DON'T assume gene birth rate = gene death rate

# Basically, we use the same scripts as before, but instead of using "lambda"
# We write lambdamu. Super easy!

cd cafe_runs/lambdamu/r1

nano newgenome_lambdamu_cafe_r1.sh

#!/mnt/griffin/handor/software/cafe
tree (Pluxyl:127.3,((Adohon:60,Cydspl:60):55.1,(((Papgla:23.8,(Pappol:14.5,(Papmac:10.8,Papxut:10.8):3.7):9.3):74.5,(Leracc:92.4,((Phosen:51.7,Pienap:51.7):33.2,(Calcec:78.6,(Danple:69.4,(Bicany:63.8,((Helmel:11.8,(Heldem:1,Hellat:1):10.8):43,((Juncoe:29.6,Melcin:29.6):8.2,Vancar:37.8):17):9):5.6):9.2):6.3):7.5):5.9):9.3,(Hypkah:105.2,((Ostnub:91.6,(Amytra:33,Ploint:33):58.6):9.7,(((Bommor:70.3,(Hylves:42.8,Mansex:42.8):27.5):18.4,Opebru:88.7):2.5,(Arcpla:68.5,(Spofru:60.6,Tricni:60.6):7.9):22.7):10.1):3.9):2.4):7.5):12.2)
load -i ../../../all_v_all_blast/output/filtered_new_genomes_Aug23_cafe_input.txt -t 20 -l ../../cafe_reports/lambdamu_r1.log -p 0.05 -filter
errormodel -model /mnt/griffin/handor/2023_CAFE_redo/cafe_runs/cafe_reports/new_genomes_errormodel_filter_r1/cafe_errormodel_0.100390625.txt -all
lambdamu -s
report ../../cafe_reports/lambdamu_r1


##### Accounting for large families ####

# If you recall, we filtered out families where the average number of copies per species was greater than 100.
# What happens if we want to put those families back in the analysis?

cd /mnt/griffin/handor/2023_CAFE_redo/cafe_runs/largefams

nano newgenome_largefam_lambdamu_cafe_r1.sh

#!/mnt/griffin/handor/software/cafe
tree (Pluxyl:127.3,((Adohon:60,Cydspl:60):55.1,(((Papgla:23.8,(Pappol:14.5,(Papmac:10.8,Papxut:10.8):3.7):9.3):74.5,(Leracc:92.4,((Phosen:51.7,Pienap:51.7):33.2,(Calcec:78.6,(Danple:69.4,(Bicany:63.8,((Helmel:11.8,(Heldem:1,Hellat:1):10.8):43,((Juncoe:29.6,Melcin:29.6):8.2,Vancar:37.8):17):9):5.6):9.2):6.3):7.5):5.9):9.3,(Hypkah:105.2,((Ostnub:91.6,(Amytra:33,Ploint:33):58.6):9.7,(((Bommor:70.3,(Hylves:42.8,Mansex:42.8):27.5):18.4,Opebru:88.7):2.5,(Arcpla:68.5,(Spofru:60.6,Tricni:60.6):7.9):22.7):10.1):3.9):2.4):7.5):12.2)
load -i /mnt/griffin/handor/2023_CAFE_redo/all_v_all_blast/output/not_sizefiltered_new_genomes_Aug23_cafe_input.txt -t 30 -l ../cafe_reports/largefams_lambdamu_r1.log -p 0.05 -filter
errormodel -model /mnt/griffin/handor/2023_CAFE_redo/cafe_runs/cafe_reports/new_genomes_errormodel_filter_r1/cafe_errormodel_0.100390625.txt -all
lambdamu -l 0.00288653 -m 0.00095399
report ../../cafe_reports/largefams_lambdamu_r1
