#!/mnt/griffin/handor/software/cafe
load -i ../all_v_all_blast/output/filtered_new_genomes_Aug23_cafe_input.txt -t 20 -l cafe_reports/butterflyrate_errormodel_building.log -filter
tree (Pluxyl:127.3,((Adohon:60,Cydspl:60):55.1,(((Papgla:23.8,(Pappol:14.5,(Papmac:10.8,Papxut:10.8):3.7):9.3):74.5,(Leracc:92.4,((Phosen:51.7,Pienap:51.7):33.2,(Calcec:78.6,(Danple:69.4,(Bicany:63.8,((Helmel:11.8,(Heldem:1,Hellat:1):10.8):43,((Juncoe:29.6,Melcin:29.6):8.2,Vancar:37.8):17):9):5.6):9.2):6.3):7.5):5.9):9.3,(Hypkah:105.2,((Ostnub:91.6,(Amytra:33,Ploint:33):58.6):9.7,(((Bommor:70.3,(Hylves:42.8,Mansex:42.8):27.5):18.4,Opebru:88.7):2.5,(Arcpla:68.5,(Spofru:60.6,Tricni:60.6):7.9):22.7):10.1):3.9):2.4):7.5):12.2)
lambda -s -t (1,((1,1)1,(((2,(2,(2,2)2)2)2,(2,((2,2)2,(2,(2,(2,((2,(2,2)2)2,((2,2)2,2)2)2)2)2)2)2)2)2,(1,((1,(1,1)1)1,(((1,(1,1)1)1,1)1,(1,(1,1)1)1)1)1)1)1)1)
report cafe_reports/butterflyrate_building_errormodel_filter_r1

# Run with
python2 /mnt/griffin/handor/software/CAFE/cafe/caferror.py -i butterflyrate_errormodel_cafe.sh -d cafe_reports/butterflyrate_errormodel_filter_r1 -v 0 -f 1

# Gave
# *Error 5: Only single lambda searches are recommended for estimating error. For more info, see the CAFE Manual.

# So let's just do a normal run with the single lambda error model

nano butterflyrate_lambdamu_cafe_r1.sh
#!/mnt/griffin/handor/software/cafe
tree (Pluxyl:127.3,((Adohon:60,Cydspl:60):55.1,(((Papgla:23.8,(Pappol:14.5,(Papmac:10.8,Papxut:10.8):3.7):9.3):74.5,(Leracc:92.4,((Phosen:51.7,Pienap:51.7):33.2,(Calcec:78.6,(Danple:69.4,(Bicany:63.8,((Helmel:11.8,(Heldem:1,Hellat:1):10.8):43,((Juncoe:29.6,Melcin:29.6):8.2,Vancar:37.8):17):9):5.6):9.2):6.3):7.5):5.9):9.3,(Hypkah:105.2,((Ostnub:91.6,(Amytra:33,Ploint:33):58.6):9.7,(((Bommor:70.3,(Hylves:42.8,Mansex:42.8):27.5):18.4,Opebru:88.7):2.5,(Arcpla:68.5,(Spofru:60.6,Tricni:60.6):7.9):22.7):10.1):3.9):2.4):7.5):12.2)
load -i ../all_v_all_blast/output/filtered_new_genomes_Aug23_cafe_input.txt -t 20 -l ./cafe_reports/butterflyrate_lambdamu_r1.log -p 0.05 -filter
errormodel -model /mnt/griffin/handor/2023_CAFE_redo/cafe_runs/cafe_reports/new_genomes_errormodel_filter_r1/cafe_errormodel_0.100390625.txt -all
lambdamu -s -t (1,((1,1)1,(((2,(2,(2,2)2)2)2,(2,((2,2)2,(2,(2,(2,((2,(2,2)2)2,((2,2)2,2)2)2)2)2)2)2)2)2,(1,((1,(1,1)1)1,(((1,(1,1)1)1,1)1,(1,(1,1)1)1)1)1)1)1)1)
report ./cafe_reports/butterflyrate_lambdamu_r1

# If we want to compare models, we have to run them with simulated data

nano compare_models.sh
#!/mnt/griffin/handor/software/cafe
tree (Pluxyl:127.3,((Adohon:60,Cydspl:60):55.1,(((Papgla:23.8,(Pappol:14.5,(Papmac:10.8,Papxut:10.8):3.7):9.3):74.5,(Leracc:92.4,((Phosen:51.7,Pienap:51.7):33.2,(Calcec:78.6,(Danple:69.4,(Bicany:63.8,((Helmel:11.8,(Heldem:1,Hellat:1):10.8):43,((Juncoe:29.6,Melcin:29.6):8.2,Vancar:37.8):17):9):5.6):9.2):6.3):7.5):5.9):9.3,(Hypkah:105.2,((Ostnub:91.6,(Amytra:33,Ploint:33):58.6):9.7,(((Bommor:70.3,(Hylves:42.8,Mansex:42.8):27.5):18.4,Opebru:88.7):2.5,(Arcpla:68.5,(Spofru:60.6,Tricni:60.6):7.9):22.7):10.1):3.9):2.4):7.5):12.2)
load -i ../all_v_all_blast/output/filtered_new_genomes_Aug23_cafe_input.txt -t 30 -l ./cafe_reports/model_comp.log -filter
lambda -l 0.00213903
genfamily butterfly_sims/ rnd -t 100
lhtest -d butterfly_sims -t (1,((1,1)1,(((2,(2,(2,2)2)2)2,(2,((2,2)2,(2,(2,(2,((2,(2,2)2)2,((2,2)2,2)2)2)2)2)2)2)2)2,(1,((1,(1,1)1)1,(((1,(1,1)1)1,1)1,(1,(1,1)1)1)1)1)1)1)1) -l 0.00213903 -o ./cafe_reports/model_comp

# A few days later, you can check the results
# This is after only 75 runs to test

# You need the log liklihood values from the two model types from your intial runs
# When calling the R script

nano lhtest.R

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
folder <- dirname(args[1])
df <- read.table(args[1])
df[,3] <- 2*(df[,1]-df[,2])
names(df) <- c('global', 'multi', 'diff')

global.lnk <- as.numeric(args[2])
multi.lnk <- as.numeric(args[3]) 

print(df)

h <- ggplot(df, aes(x=diff)) + geom_histogram() + xlim(-12, 2) + xlab("2*(global_lnL-multi_lnL)") + ylab("Count") +
    theme_bw() +
    theme(# first 5 to make background white and draw axes lines
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(),
          # now making stuff transparent
          plot.background = element_rect(fill = "transparent"),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          ## axis.ticks = element_line(),
          axis.title.x = element_text(size=12),
          axis.title.y = element_text(size=12),
          legend.background=element_rect(fill = "transparent"),
          legend.title=element_text(size=12),
          legend.text=element_text(size=12)
    )

pdf(paste0(folder, "/", "lk_null.pdf"), width=4, height=4)
h + annotate("text", -8, 10, label=paste0("p-value = \n(Counts < ", 2*(global.lnk-multi.lnk), ")/100"), size=3)
dev.off()


cut -f 2,4 ./cafe_reports/model_comp > lh_test_75run_input.txt
Rscript lhtest.R lh_test_75run_input.txt -301472.521 -296655.613

Output format =  (node ID, node ID): (0,5) (2,4) (3,37) (6,9) (8,11) (10,12) (7,15) (14,19) (16,18) (17,21) (20,23) (22,25) (24,31) (26,29) (28,30) (27,35) (32,34) (33,36) (13,39) (38,45) (40,43) (42,44) (41,53) (46,49) (48,50) (47,52) (51,55) (54,57) (56,58)
'ID'    'Newick'        'Family-wide P-value'   'Viterbi P-values'      'cut P-value'   'Likelihood Ratio'
264     (Pluxyl_9:127.3,((Adohon_5:60,Cydspl_7:60)_6:55.1,(((Papgla_5:23.8,(Pappol_4:14.5,(Papmac_5:10.8,Papxut_3:10.8)_4:3.7)_4:9.3)_4:74.5,(Leracc_4:92.4,((Phosen_5:51.7,Pienap_4:51.7)_4:33.2,(Calcec_6:78.6,(Danple_4:69.4,(Bicany_4:63.8,((Helmel_4:11.8,(Heldem_4:1,Hellat_6:1)_4:10.8)_4:43,((Juncoe_3:29.6,Melcin_5:29.6)_4:8.2,Vancar_4:37.8)_4:17)_4:9)_4:5.6)_4:9.2)_4:6.3)_4:7.5)_4:5.9)_4:9.3,(Hypkah_9:105.2,((Ostnub_8:91.6,(Amytra_7:33,Ploint_6:33)_6:58.6)_6:9.7,(((Bommor_8:70.3,(Hylves_5:42.8,Mansex_4:42.8)_5:27.5)_6:18.4,Opebru_7:88.7)_6:2.5,(Arcpla_9:68.5,(Spofru_18:60.6,Tricni_8:60.6)_8:7.9)_8:22.7)_6:10.1)_6:3.9)_6:2.4)_6:7.5)_6:12.2)_6    0.001   ((0.331414,0.604513),(0.138963,0.479976),(0.801311,0.564963),(0.22303,0.572549),(0.606327,0.526044),(0.108015,0.033047),(0.814436,0.542334),(0.84155,0.55782),(0.400493,0.763686),(0.703319,0.550174),(0.171549,0.572549),(0.805231,0.542334),(0.79298,0.572549),(0.586562,0.579643),(0.508905,5.97172e-05),(0.739728,0.624694),(0.111666,0.26813),(0.565277,0.718855),(0.00197025,0.519822),(0.271901,0.529339),(0.377926,0.809129),(0.316545,0.725243),(0.581407,0.589311),(0.298847,0.111488),(0.733881,0.0720145),(0.645686,0.863052),(0.519822,0.0322272),(0.860968,0.584063),(1.74244e-05,0.853965))
