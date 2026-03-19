library(ggsci) # if not using this, make sure to switch out the palatte in the figures
library(pairwiseAdonis)

#############################################################################################
####################################   CSS TRANSFORM    ###################################
sample_data(data_AMR.ps)$Dilution <- factor(sample_data(data_AMR.ps)$Dilution, levels = c("Standard","1:2 dilution",
                                                                                        "1:3 dilution", "1:4 dilution","1:5 dilution"))

sample_data(data_AMR.ps)$SampleGroup <- as.factor(sample_data(data_AMR.ps)$SampleGroup)


data_AMR.ps <- prune_taxa(taxa_sums(data_AMR.ps) > 0, data_AMR.ps)
any(taxa_sums(data_AMR.ps)==0) # QUADRUPLE CHECKING - nope good.

data_AMR.ps.css <- phyloseq_transform_css(data_AMR.ps, log = F)
data_AMR.ps.css.df <- as(sample_data(data_AMR.ps.css), "data.frame")



###
####
##### All samples- Beta diversity #######
####
### 
data_AMR.ps.css <- phyloseq_transform_css(data_AMR.ps, log = F)
data_AMR.ps.css.df <- as(sample_data(data_AMR.ps.css), "data.frame")


ord_explore(data_AMR.ps)

data_AMR.ps.css %>%
  tax_transform(rank = "unique", trans = "hellinger") %>%
  dist_calc(dist = "bray") %>%
  ord_calc(
    method = "auto"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "SampleGroup", fill = "SampleGroup",
    shape = "LIbraryType", alpha = 0.5,
    size = 2
  ) + 
  theme_classic() +
  theme(plot.caption = element_blank()) +
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse(
    ggplot2::aes(colour = SampleGroup)
  )


data_AMR.ps.css %>%
  tax_transform(rank = "unique", trans = "hellinger") %>%
  dist_calc(dist = "bray") %>%
  ord_calc(
    method = "auto"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Dilution", fill = "Dilution",
    shape = "LIbraryType", alpha = 0.5,
    size = 2
  ) + 
  theme_classic() +
  theme(plot.caption = element_blank()) +
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse(
    ggplot2::aes(colour = Dilution)
  )

#### Ordination storage

data_AMR.ord <- data_AMR.ps.css %>%
  tax_transform(rank = "unique", trans = "hellinger") %>%
  dist_calc(dist = "bray") %>%
  ord_calc(
    method = "auto"
  )

data_AMR.dist <- data_AMR.ps.css %>%
  tax_transform(rank = "unique", trans = "hellinger") %>%
  dist_calc(dist = "bray") 

data_AMR.dist@dist

## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
SampleGroup.permanova <- pairwise.adonis(data_AMR.dist@dist, data_AMR.ps.css.df$SampleGroup, perm = 9999, p.adjust.m = "BH")
SampleGroup.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
SampleGroup.disper <- betadisper(data_AMR.dist@dist, data_AMR.ps.css.df$SampleGroup)
SampleGroup.permdisp <- permutest(SampleGroup.disper, permutations = 9999, pairwise = T)
SampleGroup.permdisp # looks like a few are significant

## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
OriginalSample.permanova <- pairwise.adonis(data_AMR.dist@dist, data_AMR.ps.css.df$OriginalSample, perm = 9999, p.adjust.m = "BH")
OriginalSample.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
OriginalSample.disper <- betadisper(data_AMR.dist@dist, data_AMR.ps.css.df$OriginalSample)
OriginalSample.permdisp <- permutest(OriginalSample.disper, permutations = 9999, pairwise = T)
OriginalSample.permdisp # looks like a few are significant


## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
Dilution.permanova <- pairwise.adonis(data_AMR.dist@dist, data_AMR.ps.css.df$Dilution, perm = 9999, p.adjust.m = "BH")
Dilution.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
Dilution.disper <- betadisper(data_AMR.dist@dist, data_AMR.ps.css.df$Dilution)
Dilution.permdisp <- permutest(Dilution.disper, permutations = 9999, pairwise = T)
Dilution.permdisp # looks like a few are significant

## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
InputDNA.permanova <- pairwise.adonis(data_AMR.dist@dist, data_AMR.ps.css.df$InputDNA, perm = 9999, p.adjust.m = "BH")
InputDNA.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
InputDNA.disper <- betadisper(data_AMR.dist@dist, data_AMR.ps.css.df$InputDNA)
InputDNA.permdisp <- permutest(InputDNA.disper, permutations = 9999, pairwise = T)
InputDNA.permdisp # looks like a few are significant





### Now with CSS values ####
data_AMR.ps.css

# Bray curtis distance matrix did converge with all samples
# Had to use a different beta diversity index, euclidean, calculated on hellinger transformed counts
data.dist <- vegdist(decostand(t(otu_table(data_AMR.ps.css)), "hell"), "euclidean") 
data.ord <- vegan::metaMDS(comm = t(data.dist), distance = "none", try = 10, trymax = 999, autotransform = F)

## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
SampleGroup.permanova <- pairwise.adonis(data_AMR.dist@dist, data_AMR.ps.css.df$SampleGroup, perm = 9999, p.adjust.m = "BH")
SampleGroup.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
SampleGroup.disper <- betadisper(data_AMR.dist@dist, data_AMR.ps.css.df$SampleGroup)
SampleGroup.permdisp <- permutest(SampleGroup.disper, permutations = 9999, pairwise = T)
SampleGroup.permdisp # looks like a few are significant


###
####
##### All samples- Beta diversity #######
####
### 

# Bray curtis distance matrix did converge with all samples
# Had to use a different beta diversity index, euclidean, calculated on hellinger transformed counts
data.dist <- vegdist(decostand(t(otu_table(ps_class.css)), "hell"), "euclidean") 
data.ord <- vegan::metaMDS(comm = t(data.dist), distance = "none", try = 10, trymax = 999, autotransform = F)
#plot_ordination(data.css, data.ord, color = "Dilution") 
plt_ord_by_Dilution <- plot_ordination(ps_class.css, data.ord, color = "Dilution") +
  theme_bw() +
  labs(title ="Sample type") +
  stat_ellipse(aes(fill= Dilution), geom="polygon", alpha = 0.25) +
  theme(legend.position = "right",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", size = 1.0),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
plt_ord_by_Dilution

# permanova test of beta diversity, by Dilution
data.adonis <- adonis2(data.dist ~ Dilution, data = data_noSNP.css.df, permutations = 9999)
data.adonis # 7e-04, R2 =0.13205 

## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
combined_order.permanova <- pairwise.adonis(data.dist, data_noSNP.css.df$Dilution, perm = 9999, p.adjust.m = "BH")
combined_order.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
combined_order.disper <- betadisper(data.dist, data_noSNP.css.df$Dilution)
Dilution.permdisp <- permutest(combined_order.disper, permutations = 9999, pairwise = T)
Dilution.permdisp # looks like a few are significant


# This requires a directory called figures to be in your working directory
# All 3 lines have to be run for the file to be created.
png("Figures/Allsamples_NMDS_Dilution.png", width = 1200, height = 800)
plt_ord_by_Dilution
dev.off()

# svg format is nice because it stores all graphic data and can be requested by some manuscripts
svg("Figures/Allsamples_NMDS_Dilution.svg", width = 500, height = 500)
plt_ord_by_Dilution
dev.off()

## STATS BETWEEN ALL GROUPS
# create distance matrix we use to ordinate (Bray-Curtis)
noSNP.dist <- vegdist(t(otu_table(data_AMR.ps.css)), method = "bray")

## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
all.groups.permanova <- pairwise.adonis(noSNP.dist, data_AMR.ps.css.df$Dilution, perm = 9999, p.adjust.m = "BH")
all.groups.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
all.groups.disper <- betadisper(noSNP.dist, data_AMR.ps.css.df$Dilution)
all.groups.permdisp <- permutest(all.groups.disper, permutations = 9999, pairwise = T)
all.groups.permdisp # looks like a few are significant
