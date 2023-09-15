# AGAT iso dataset Bayesian model, with corrected host order information
# September 2023 version

## Prep ##
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(tidyverse)
library(brms)
library(modelr)
library(tidybayes)
library(data.table)
library(cmdstanr)
library(dplyr)

theme_set(theme_minimal())

setwd("/mnt/griffin/handor/2023_CAFE_redo/")

###


d <- data.table::fread('/mnt/griffin/handor/2023_CAFE_redo/cafe_runs/cafe_reports/R_ready_orders_appended_counts_new_genomes_lambdamu_p05_r1.tsv') %>%
  mutate(
    Species = str_trim(Species),
    phylo = Species,
    # this line below  will need to be changed depending on input
    across(Fam_11:Fam_22837)
  ) %>%
  pivot_longer(starts_with('Fam'), names_to = 'family', values_to = 'num_genes') %>%                                                                 
  rename_with(tolower)

###
# Load the tree:
tree <- ape::read.tree('tree/cydspl_inc_tree.tre')

# We don't have host order data for Hypkah, so we drop it here
tree <- ape::drop.tip(tree, 'Hypkah')
tree_cov <- ape::vcv(tree)

# So, there is actually a large number of gene families that have very few genes in them, and not
# much variation between species. For very small families, it won't really be feasible to detect
# a relationship with host range.

pdf("mean_fam_size_plot.pdf")
d %>% group_by(family) %>% summarise(variance = var(num_genes)) %>%
  ggplot(aes(variance)) + geom_histogram(bins = 200) + xlim(0, 15)

d %>% group_by(family) %>% summarise(total = mean(num_genes)) %>%
  ggplot(aes(total)) + geom_histogram(bins = 200) + xlim(0, 15) + geom_vline(xintercept = 2, lty = 2)

dev.off()

# So I propose we fit our model on a subset of the data, only including reasonably large families.
# This will make the model much easier to fit, because there is just much less data, but also
# because the fitting problems were worst at low gene counts. It's difficult to say what a good
# cut-off is; here, we pick 2 or more genes per species. You can change this to a
# different cut-off. I do recommend filtering on the average number of genes, not on the variance.

d_2gene <- d %>% group_by(family) %>% filter(mean(num_genes) >= 2)

set_cmdstan_path(path="/home/handor/.cmdstan/cmdstan-2.32.0")

m_2gene <- brm(
  ## Defining our model:
  # formula:
  num_genes ~ 1 + num_orders + (1 + num_orders | family) +
    (1 | species) + (1 | gr(phylo, cov = A)),
  # data:
  data = d_2gene, data2 = list(A = tree_cov),
  family = 'negbinomial',
  ## Our sampling strategy across chains:
  warmup = 1000, iter = 3000, # so 2000 actual samples per chain
  # I think 4 chains is good here, because there is a lot of variation in sampling speed between
  # two chains. The 8000 samples I got here was more than enough for good effective sample sizes
  # (ESS), and Rhat values very close to 1. Rhat compares the results from the different chains to
  # see if we have good convergence, if it's 1 then the chains are giving the same results.
  # You can speed up fitting by increasing the number of threads here. Just note that this is the
  # number of threads _per chain_, so here it will use 4 * 6 = 24 threads. Speed increase is not
  # linear though.
  cores = 4, chains = 4, threads = 6,
  ## Other options
  # using cmdstanr is often much faster (sometimes 3x times!) and is needed to use threads.
  backend = 'cmdstanr',
  # I got a warning about divergent transititions. Increasing adapt_delta to be closer to 1 will
  # help, at the cost of slower sampling (default is 0.8).
  control = list(adapt_delta = 0.99, max_treedepth = 15), save_model ='Sep2023_AGATiso_2gene_reduced_set_model_nb.rds' 
)

# My timings for the fit:
# Mean chain execution time: 4073.5 seconds.
# Total execution time: 4700.8 seconds.


# look at the posterior predictive checks:
pdf("ppcheck.pdf")
pp_check(m_2gene, ns = 50) + xlim(-10, 125)
pp_check(m_2gene, ns = 50) + xlim(-1, 10)
dev.off()

# So still not ideal, but I think this is just about good enough.

# Full model summary:
summary(m_2gene)

# save it
saveRDS(m_2gene, file="m_2gene_AGATiso_newgenomes.RData")
# does it load? Yes.
test<-readRDS("m_2gene_AGATiso_newgenomes.RData")

# So now that we have a model, let's get the per-family estimates of the slope between the gene
# family size and the host range:
pd2 <- spread_draws(m_2gene, b_num_orders, r_family[family, parameter]) %>%
  filter(parameter == 'num_orders') %>%
  mutate(slope = b_num_orders + r_family)

# And plot:
plot1_99 <- pd2 %>%
  mutate(incidence_ratio = -((1 - exp(slope)) * 100)) %>%
  ggplot(aes(incidence_ratio, fct_reorder(family, incidence_ratio))) +
  # You have the change the size parameter here as you add more
  stat_interval(size = 0.5) +
  geom_vline(xintercept = 0, lty = 1) +
  scale_color_brewer(labels = function(x) scales::percent(as.numeric(x), accuracy = 1)) +
  labs(
    x = 'Effect of host range on gene family size\n(% change for 1 additional host order)',
    y = 'Gene family',
    color = 'Credible Interval:'
  ) +
  theme_minimal() +
  theme(legend.position = 'top', panel.grid.major.y = element_blank(), axis.text.y = element_blank())
ggsave('newgenomes_plot1_99_nogray.png', plot1_99, width = 5, height = 10)

# So as you see, this is a much smaller set of families. We do see a number of them where we can
# clearly see that the slope is not zero. So let's get a list of those families, and plot the model
# fits for those family (with the raw data) to get an idea whether our model is detecting things
# that indeed make sense.

# Summarise the posteriors to get estimates and credible intervals:
family_summ <- median_qi(pd2)

# Filter to get only the families with significant slope:
the_chosen <- filter(family_summ, sign(slope.lower) == sign(slope.upper))

# Now get the full regression lines with uncertainty:
nd <- expand_grid(
  family = the_chosen$family,
  num_orders = seq(1, 6, l = 60)
)
nd$num_genes <- m_2gene %>%
  posterior_epred(newdata = nd, re_formula = ~ (1 + num_orders | family)) %>%
  as.data.frame() %>% as.list()

nd <- unnest(nd, num_genes) %>%
  group_by(family, num_orders) %>%
  median_qi()

family_plot <- function(raw_data, regression_lines) {
  ggplot(mapping = aes(num_orders, num_genes)) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), regression_lines, alpha = 0.3) +
    geom_line(data = regression_lines, color = 'firebrick') +
    geom_count(data = raw_data, alpha = 0.6) +
    scale_size_area(max_size = 3) +
    scale_x_continuous(breaks = 1:6) +
    facet_wrap(~family, scales = 'free') +
    labs(
      x = 'Number of host orders',
      y = 'Gene family size',
      size = 'Count'
    ) +
    theme_bw()
}
family_plot

# Plot the families with positive slopes
plot2_99 <- family_plot(
  filter(d_2gene, family %in% filter(the_chosen, slope > 0)$family),
  filter(nd, family %in% filter(the_chosen, slope > 0)$family)
)
ggsave('newgenome_AGATiso_plot2_99.png', plot2_99, width = 7, height = 5)

# Plot the families with negative slopes
plot3_99 <- family_plot(
  filter(d_2gene, family %in% filter(the_chosen, slope <= 0)$family),
  filter(nd, family %in% filter(the_chosen, slope <= 0)$family)
)
ggsave('newgenome_AGATiso_plot3_99.png', plot3_99, width = 7, height = 6)

# redo this, but with annotation labels

# Filter out families that are unknown/TEs
messyfams<-c("Fam_358","Fam_1031","Fam_42","Fam_52","Fam_99","Fam_105","Fam_116","Fam_186","Fam_301")
`%notin%` <- Negate(`%in%`)
filtered_chosen<-the_chosen[the_chosen$family %notin% messyfams,]

plot2_99 <- family_plot(
  filter(d_2gene, family %in% filter(filtered_chosen, slope > 0)$family),
  filter(nd, family %in% filter(filtered_chosen, slope > 0)$family)
)

genefam_labs_plot2<- c("cuticle protein", "K02A2.6-like", "serine-type \nendopeptidase activity","Alcohol dehydrogenase \ntranscription factor \nand/or MADF","ubiquitin-dependent\nendocytosis")
names(genefam_labs_plot2) <- c("Fam_138","Fam_161","Fam_315","Fam_492","Fam_527")

color_transparent1 <- adjustcolor("#44AA99", alpha.f = 0.50)

plot2_99_relabel <- plot2_99 + facet_wrap(~family, ncol = 3, scales = "free_y", labeller = labeller(family = genefam_labs_plot2)) + theme(strip.text= element_text(size =9), strip.background = element_rect(fill = color_transparent1))

ggsave('newgenomes_plot2_99_relabel.png', plot2_99_relabel, width = 8, height = 8)

###
plot3_99 <- family_plot(
  filter(d_2gene, family %in% filter(filtered_chosen, slope <= 0)$family),
  filter(nd, family %in% filter(filtered_chosen, slope <= 0)$family) 
)
genefam_labs_plot3<- c("chorion protein", "chorion protein", "trypsin-like \nserine protease", "trypsin-like \nserine protease")
names(genefam_labs_plot3) <- c("Fam_76", "Fam_123", "Fam_140", "Fam_204")

color_transparent2 <- adjustcolor("#44AA99", alpha.f = 0.50)

plot3_99_relabel <- plot3_99 + facet_wrap(~family, ncol = 3, nrow=2, scales = "free_y", labeller = labeller(family = genefam_labs_plot3)) + theme(strip.text= element_text(size =9), strip.background = element_rect(fill = color_transparent2))

ggsave('newgenomes_plot3_99_relabel.png', plot3_99_relabel, width = 8, height = 8)

# Make it prettier
stacks <- ggarrange (NULL, plot2_99_relabel, NULL, plot3_99_relabel, ncol=1, heights =c(0.1, 1.25,0.1,1.25), widths=1.5, common.legend=TRUE, labels =c("4b: increasing with generalization", NA,"4c: increasing with specialization", NA),vjust=0.55,hjust=0, legend = "bottom")
stacks2 <- ggarrange(NULL, plot1_99, nrow = 2, ncol = 1, heights = c(0.05, 2.05), widths = 1, labels = c("4a", NA),vjust=0.8)
attempt3<-ggarrange(NULL, NULL, stacks2, stacks, nrow = 2, ncol = 2, heights = c(0.05, 1), widths = c(1,2))

pdf(file="newgenomes_combined_sigfams2.pdf", width = 12, height = 8)
attempt3
dev.off()
# ggsave('corrected_AGAT_iso_combined_sigfams2.png', attempt3, width = 10, height = 8)

#######################################

# What about the dataset that includes gene families > 100 copies per species?
large_data <- data.table::fread('/mnt/griffin/handor/2023_CAFE_redo/cafe_runs/cafe_reports/R_ready_orders_appended_counts_largefams_lambdamu_p05_r1.tsv') %>%
  mutate(
    Species = str_trim(Species),
    phylo = Species,
    # This line below  will need to be changed depending on input
    across(Fam_1:Fam_22837)
  ) %>%
  pivot_longer(starts_with('Fam'), names_to = 'family', values_to = 'num_genes') %>%                                                                 
  rename_with(tolower)

large_data %>% group_by(family) %>% summarise(variance = var(num_genes)) %>%
  ggplot(aes(variance)) + geom_histogram(bins = 200) + xlim(0, 15)
large_data %>% group_by(family) %>% summarise(total = mean(num_genes)) %>%
  ggplot(aes(total)) + geom_histogram(bins = 200) + xlim(0, 15) + geom_vline(xintercept = 2, lty = 2)

# Subset down to families where average gene count two or more.
large_filtered_fams <- large_data %>% group_by(family) %>% filter(mean(num_genes) >= 2)

set_cmdstan_path(path="/home/handor/.cmdstan/cmdstan-2.32.0")

large_fams_model_sizefilter <- brm(
  ## Defining our model:
  # formula:
  num_genes ~ 1 + num_orders + (1 + num_orders | family) +
    (1 | species) + (1 | gr(phylo, cov = A)),
  # data:
  data = large_filtered_fams, data2 = list(A = tree_cov),
  family = 'negbinomial',
  ## Our sampling strategy across chains:
  warmup = 1000, iter = 3000, # so 2000 actual samples per chain
  cores = 4, chains = 4, threads = 6,
  backend = 'cmdstanr',
  control = list(adapt_delta = 0.99, max_treedepth = 15), save_model ='largefams_2gene_reduced_set_model_nb.rds' 
)

saveRDS(large_fams_model_sizefilter, file="m_largefam_2gene_AGATiso_newgenomes.RData")

# I cleared my environment here to prevent accidentally reading in any weird files.
# Meant I needed to reload in a few dfs above though.

m_largefam <- readRDS("m_largefam_2gene_AGATiso_newgenomes.RData")

lf2 <- spread_draws(m_largefam, b_num_orders, r_family[family, parameter]) %>%
  filter(parameter == 'num_orders') %>%
  mutate(slope = b_num_orders + r_family)

large_family_summ <- median_qi(lf2)

# Filter to get only the families with significant slope:
lf_chosen <- filter(large_family_summ, sign(slope.lower) == sign(slope.upper))

# Now get the full regression lines with uncertainty:
lf_nd <- expand_grid(
  family = lf_chosen$family,
  num_orders = seq(1, 6, l = 60)
)

lf_nd$num_genes <- m_largefam %>%
  posterior_epred(newdata = lf_nd, re_formula = ~ (1 + num_orders | family)) %>%
  as.data.frame() %>% as.list()

nd_unnest <- unnest(lf_nd, num_genes) %>%
  group_by(family, num_orders) %>%
  median_qi()

lf_family_plot <- function(raw_data, regression_lines) {
  ggplot(mapping = aes(num_orders, num_genes)) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), regression_lines, alpha = 0.3) +
    geom_line(data = regression_lines, color = 'firebrick') +
    geom_count(data = raw_data, alpha = 0.6) +
    scale_size_area(max_size = 3) +
    scale_x_continuous(breaks = 1:6) +
    facet_wrap(~family, scales = 'free') +
    labs(
      x = 'Number of host orders',
      y = 'Gene family size',
      size = 'Count'
    ) +
    theme_bw()
}

# Plot the families with positive slopes
lf_plot2_99 <- lf_family_plot(
  filter(large_filtered_fams, family %in% filter(lf_chosen, slope > 0)$family),
  filter(nd_unnest, family %in% filter(lf_chosen, slope > 0)$family)
)
ggsave('lf_newgenome_AGATiso_plot2_99.png', lf_plot2_99, width = 7, height = 5)

# Plot the families with negative slopes
lf_plot3_99 <- lf_family_plot(
  filter(large_filtered_fams, family %in% filter(lf_chosen, slope <= 0)$family),
  filter(nd_unnest, family %in% filter(lf_chosen, slope <= 0)$family)
)
ggsave('lf_newgenome_AGATiso_plot3_99.png', lf_plot3_99, width = 7, height = 6)