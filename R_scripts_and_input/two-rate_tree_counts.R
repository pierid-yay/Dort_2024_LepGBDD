# Feb 2024 #

# Regenerating figures, in light of two-rate models

# First, make sure everything is nice and cleared.
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

# My installation of R on the server is not cooperating
# So let's reinstall the needed packages on my local computer.

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("treeio")
# BiocManager::install("ggtree")
# install.packages("phytools")

# Load in our libraries
library(ggplot2)
library(ape)
library(phytools)
library(ggrepel)
library(ggpubr)
library(scales)
library(ggtree)
library(ggstance)

# Load in / set the three data frames we need and set the x max variable:
setwd ("C:/Users/hanna/Desktop/Ch1_tree_redo")

# The file with manually set ggtree nodes
## Note: Before running any analyses, a column must be added to the .csv file with the ggtree
## node labels. Unfortunately this is done manually by comparing a figure with the ggtree
## labels to a figure with the CAFE labels (visualized with something like seaview or Figtree).
## If this isn't done, data will map to the wrong branch in the tree.
dataframe <- read.csv("feb_newgenome_tree_nodes_and_counts.csv", sep=",")

# The tree file
tree_file = "cydspl_inc_tree.tre"
tree = read.tree(tree_file)
node_check = T

# And the file with information on the number of host orders
tips_info=read.csv("new_genome_hostcounts_forR.csv")
head(tips_info)

# Scale xmax depending on the length of your tree.
xmax = 200

#########################
# Node label checking -- this block displays a figure with the ggtree node labels.
if(node_check){
  node_test = ggtree(tree, size=1, ladderize=F) +
    ggplot2::xlim(0, xmax) +
    geom_tiplab(color="#333333", fontface='italic', size=5) +
    geom_text(aes(label=node), hjust=-.3, vjust=-.3, color="#ff6db6") +
    geom_nodepoint(color="#666666", alpha=0.85, size=4)
  print(node_test)
  stop("Node check OK.")
}
################
# IMPORTANT: Once the ggtree.node column has been filled in on the .csv file,
# (which it has)
# sort the nodes
dataframe = dataframe[order(dataframe$ggtree.node),]
tips = subset(dataframe, Node.type == "Tip")
internals = subset(dataframe, Node.type == "Internal")
##########

# If plotting Fam_489 counts
# Which are the 7tm odorant receptors
# That contracted at the butterfly origin node.

tree_fig_489 = ggtree(tree, size=1, ladderize=F, aes(color=dataframe$Fam_489_count)) +
  scale_color_continuous(name='Genes in odorant \n receptor family', low="#bdd7e7", high="#3182bd") +
  ggplot2::xlim(0, xmax) +
  geom_tiplab(aes(label=dataframe$Species), color="#333333", nudge_x =1, fontface='italic', size=4) +
  geom_label(aes(x = branch, label=paste0(dataframe$Fam_489_count)), colour='black', nudge_y=0.15, size =3, label.padding = unit(0.1, "lines")) +
  theme(legend.position="bottom",
        plot.caption=element_text(hjust=0))

pdf(file ="resize_family_489_counts.pdf", height = 8, width = 8)
print(tree_fig_489)
dev.off()

# # If plotting Fam_264 counts, the OBPs
# tree_fig_264 = ggtree(tree, size=1, ladderize=F, aes(color=dataframe$Fam_264_count)) +
#   scale_color_continuous(name='Genes in odorant n/binding protein family', low="#648FFF", high="#E700B5") +
#   ggplot2::xlim(0, xmax) +
#   geom_tiplab(aes(label=dataframe$Species), color="#333333", fontface='italic', size=4) +
#   geom_label(aes(x = branch, label=paste0(dataframe$Fam_264_count)), colour='#175D68', nudge_y=0.15, size =3, label.padding = unit(0.1, "lines")) +
#   theme(legend.position="right",
#         plot.caption=element_text(hjust=0))

# pdf(file ="resize_family_264_counts.pdf", height = 8, width = 8)
# print(tree_fig_264)
# dev.off()


# I don't like how these turned out. 
# Is there a way to display both counts at once?

# Looks like we can do it, just have the tiplab only on for one tree

p2 = ggtree(tree, size=1, ladderize=F, aes(color=dataframe$Fam_264_count)) +
  scale_color_continuous(name='Genes in odorant \nbinding protein family', low="#fde0dd", high="#E700B5") +
  ggplot2::xlim(0, xmax) +
  scale_x_reverse() +
  geom_label(aes(x = branch, label=paste0(dataframe$Fam_264_count)), colour='black', nudge_y=0.15, size =3, label.padding = unit(0.1, "lines")) +
  theme(legend.position="bottom",
        plot.caption=element_text(hjust=0))

#combo <- cowplot::plot_grid(tree_fig_489, p2, ncol=2, widt)
combo2 <- ggarrange(tree_fig_489, p2, ncol=2, nrow=1, widths = c(1.5,1))
pdf(file ="combo_counts.pdf", height = 9, width = 9)
print(combo2)
dev.off()

####### Tree with species host plant order counts ######

t <-ggtree(tree, size=1, ladderize=F) +
  geom_tiplab(aes(label=dataframe$Node.ID), size=4, color="#333333", align=TRUE)

t + 
  geom_facet(panel = "Number of host plant orders", data = tips_info, geom = ggstance::geom_barh, 
             aes(x =  host_orders, color =  host_orders, fill =  host_orders), 
             stat = "identity", width = .6) +
  theme_tree2(legend.position="none")

par(mfrow=c(1,1))
pdf(file="Species_tree_w_hostorder_barplot.pdf", width=30,height=15,pointsize = 12)
t + 
  geom_facet(panel = "Number of host plant orders", data = tips_info, geom = ggstance::geom_barh, 
             aes(x =  host_orders, fill =  host_orders, color = host_orders ), 
             stat = "identity", width = .6) +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(3, "BuPu")) +
  scale_colour_gradientn(colors = RColorBrewer::brewer.pal(3, "BuPu")) +
  theme_tree2(legend.position="none")

dev.off()

# Full species names can be added in inkscape by directly modifying the pdf file
# It's easier to change things that way instead of trial-and-erroring it with font sizes

#### Tree with all expansion and contraction counts #####
# rapids in parentheses

pdf(height=13,width=10,file="new_genome_isoAGAT_lambmu.pdf")

treefig_n = ggtree(tree, size=1, ladderize=F, color = "#000000") +
  theme_tree2() +
  ggplot2::xlim(0, xmax) +
  geom_tiplab(aes(label=dataframe$Node.ID), size=3, color="#333333", fontface='italic', offset=20) +
  ggplot2::geom_label(aes(x = branch, label=paste0(dataframe$Expansions," (",dataframe$Rapid.expansions,")")), colour='#175D68', nudge_y=0.15, size =3, label.padding = unit(0.1, "lines")) +
  ggplot2::geom_label(aes(x = branch, label=paste0(dataframe$Contractions," (",dataframe$Rapid.contractions,")")), colour='#E05F72', nudge_y=-0.15, size =3, label.padding = unit(0.1, "lines")) +
  theme(legend.position="right",
        plot.caption=element_text(hjust=2)) 

treefig_n

dev.off()