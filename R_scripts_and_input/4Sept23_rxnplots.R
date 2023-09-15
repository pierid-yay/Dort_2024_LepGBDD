### Reaction norm plots and t-test ###

# I didn't like putting protein counts against the tree
# This is a solution to that, now with updated protein counts. 

# set working directory 
setwd("C:/Users/Hanna Dort/Desktop/Rscripts")

# Load in expanded counts dataframe
df <- read.csv("4Sept23_updated_long_counts_df.csv", sep=",")
head (df)
# species   type filtered count  X
# 1 Pluxyl  native   unfilt 22256 NA
# 2 Pluxyl  braker   unfilt 17590 NA
# 3 Pluxyl  braker     filt 17556 NA
# 4 Pluxyl  native     filt 13961 NA
# 5  Adohon braker   unfilt 25283 NA
# 6  Adohon braker     filt 18291 NA

library(tidyverse)
library(ggrepel)
library(dplyr)
library(ggpubr)
library(car)
library(ggstance)
library(wesanderson)


# First, we want to see if filtered braker2 dataset routinely had more proteins than filtered native dataset
# We should assess that the difference btween the two counts is normally distributed before doing this.
# With a Shapiro-Wilk test?

# Mutate the dataframe so we can see differences between filtered braker counts

filt_counts <- df %>% 
  pivot_wider(names_from = type, values_from = count) %>% 
  mutate(diff=braker-native) %>% 
  dplyr::filter(filtered=="filt")

# qqplot it for a visual check, using car package
qqPlot(filt_counts$diff)
# Looks pretty normal.

# Now for the actual test:
shapiro.test(filt_counts$diff)
# W = 0.95767, p-value = 0.4984
# So no deviation from normality!

# Paired t-test time.
res<-t.test(filt_counts$braker, filt_counts$native, paired = TRUE)
res

# data:  filt_counts$braker and filt_counts$native
# t = 10.188, df = 19, p-value = 3.894e-09
# alternative hypothesis: true mean difference is not equal to 0
# 95 percent confidence interval:
# 4450.151 6751.349
# sample estimates:
# mean difference 
# 5600.75 

### Time to get plotting! ###
df_natives<-subset(df, type=="native")
df_natives

# Subset so it's only braker2 protein set data
df_braker<-subset(df, type=="braker")
df_braker

# subset main df so it's only filtered protein set data
df_filt<-subset(df, filtered=="filt")
df_filt

# Calculate percent change caused by filtering
# For each annotation type.
b2_pivot_df <- df %>% 
  pivot_wider(names_from = filtered, values_from = count) %>% 
  mutate(b2_pct_decrease=(abs((filt-unfilt)/unfilt)*100)) %>% 
  dplyr::filter(type=="braker")

nat_pivot_df <- df %>% 
  pivot_wider(names_from = filtered, values_from = count) %>% 
  mutate(native_pct_decrease=(abs((filt-unfilt)/unfilt)*100)) %>% 
  dplyr::filter(type=="native")

testdf <- nat_pivot_df %>% select(species, native_pct_decrease)

# Place both pct differences in same df.
join_diffs <- left_join(b2_pivot_df, testdf, by='species')

# Make a palette for later
col_set <- c("#1b9e77", "#d95f02", "#7570b3", "#FF0000", "#F98400", "#5BBC6D", "#F71D7F", "#5BBCD6" )

# filtered_vs_unfiltered_native
jitter_mag <- 0.03
gg <- ggplot(data=df_natives, aes(x=filtered, y=count, group=species, col=species)) +
  geom_point() +
  geom_line(linewidth=0.5,
            alpha=0.5,
            position=position_jitter(w=0, h=jitter_mag)) +
  ylab("Protein count") +
  xlab("Filtering status of gff") +
  ylim(10000,31000)+
  ggtitle("Native annotation protein counts")+
  geom_label_repel(aes(label = ifelse(filtered=="filt", as.character(species),'')),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.linetype = "dotted",
                   nudge_x = 0.4,
                   size = 3,
                   hjust = 0,
                   direction = "y") +
  scale_color_manual(values = rep(c(col_set),10), na.value="gray") +
  guides(col = "none") +
  theme_pubr() +
  scale_x_discrete(limits=rev)+
  NULL

gg

# filtered_vs_unfiltered_braker
jitter_mag <- 0.03
gg2 <- ggplot(data=df_braker, aes(x=filtered, y=count, group=species, col=species)) +
  geom_point() +
  geom_line(linewidth=0.5,
            alpha=0.5,
            position=position_jitter(w=0, h=jitter_mag)) +
  ylab("Protein count") +
  xlab("Filtering status of gff") +
  ylim(10000,31000)+
  ggtitle("Braker2 gff protein counts")+
  geom_label_repel(aes(label = ifelse(filtered=="filt", as.character(species),'')),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.linetype = "dotted",
                   nudge_x = 0.4,
                   size = 3,
                   hjust = 0,
                   direction = "y") +
  scale_color_manual(values = rep(c(col_set),10), na.value="gray") +
  guides(col = "none") +
  theme_pubr() +
  scale_x_discrete(limits=rev)+
  NULL
gg2

# filtered b2 vs filt native
jitter_mag <- 0.03
col_set <- c("#1b9e77", "#d95f02", "#7570b3", "#FF0000", "#F98400", "#5BBC6D", "#F71D7F", "#5BBCD6" )
gg3 <- ggplot(data=df_filt, aes(x=type, y=count, group=species, col=species)) +
  geom_point() +
  geom_line(size=0.5,
            alpha=0.5,, 
            position=position_jitter(w=0, h=jitter_mag)) +
  ylab("Protein count") +
  xlab("Source of gff") +
  ylim(10000,25000)+
  scale_color_manual(values = rep(c(col_set),10), na.value="gray") +
  guides(col = "none") +
  ggtitle("Filtered gff protein counts by source")+
  geom_label_repel(aes(label = ifelse(type=="native", as.character(species),'')),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.linetype = "dotted",
                   nudge_x = 0.4,
                   size = 5,
                   hjust = 0,
                   direction = "y") +
  theme_pubr() +
  NULL
gg3

# Finally, comparing percent decrease between native and B2

# The plot
jitter_mag <- 0.03
gg4 <- ggplot(data=join_diffs, aes(x=native_pct_decrease, y=b2_pct_decrease, group=species)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "twodash", color = "darkgreen") +
  xlim(0,60)+
  ylim(0,60)+
  ylab("B2 annotation % decrease") +
  xlab("Native annotation % decrease") +
  ggtitle("Comparison of % decrease in protein counts following filtering")+
  # geom_label_repel(aes(label = species),
  #                  box.padding   = 0.35, 
  #                  point.padding = 0.5,
  #                  segment.linetype = "dotted",
  #                  nudge_x = 0.4,
  #                  size = 3,
  #                  hjust = 0) +
  scale_color_manual(values = rep(c(col_set),10), na.value="gray") +
  geom_text_repel(aes(label=species),hjust=0,vjust=0,nudge_x = 0,nudge_y = 0, segment.linetype="dotted", segment.color = "hotpink",size=3)+
  theme_pubr() +
  NULL
gg4


blep <- ggarrange(gg,gg2,gg3,gg4,ncol=2,nrow=2)
# Honestly, I don't think gg4 is that informative, 
# it should probably go in a supplement somewhere

# We need the space for the more interesting panels.
bloop <- ggarrange(gg,gg2,ncol=2, nrow=1)
blarp <- ggarrange(bloop,gg3,ncol=1, nrow=2)

pdf(file="ggarranged_rx_norms2.pdf", width=15,height=15)
blarp
dev.off()

# gg4 as a solo pdf
pdf("gg4.pdf")
gg4
dev.off()
