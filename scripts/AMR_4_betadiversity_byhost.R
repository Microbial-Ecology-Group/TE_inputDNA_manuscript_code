# AMR beta diversity by host
library(ggsci) # if not using this, make sure to switch out the palatte in the figures
library(pairwiseAdonis)

#############################################################################################
##  CSS TRANSFORM    ###################################

sample_data(group.ps)$SampleGroup <- as.factor(sample_data(group.ps)$SampleGroup)


group.ps <- prune_taxa(taxa_sums(group.ps) > 0, group.ps)
any(taxa_sums(group.ps)==0) # QUADRUPLE CHECKING - nope good.

group.ps.css <- phyloseq_transform_css(group.ps, log = F)
group.ps.css.df <- as(sample_data(group.ps.css), "data.frame")

###
####
## Split phyloseq objects by species #######
####
### 

# DatasetA 
DatasetA.group.ps.css <- subset_samples(group.ps.css, SampleGroup != "Canine fecal")
DatasetA.group.ps.css <- prune_taxa(taxa_sums(DatasetA.group.ps.css) > 0, DatasetA.group.ps.css)
any(taxa_sums(DatasetA.group.ps.css)==0) # QUADRUPLE CHECKING - nope good.
DatasetA.group.ps.css.df <- as(sample_data(DatasetA.group.ps.css), "data.frame")

# Bovine 
Bovine.group.ps.css <- subset_samples(group.ps.css, SampleGroup == "Bovine fecal")
Bovine.group.ps.css <- prune_taxa(taxa_sums(Bovine.group.ps.css) > 0, Bovine.group.ps.css)
any(taxa_sums(Bovine.group.ps.css)==0) # QUADRUPLE CHECKING - nope good.
Bovine.group.ps.css.df <- as(sample_data(Bovine.group.ps.css), "data.frame")

# Avian 
Avian.group.ps.css <- subset_samples(group.ps.css, SampleGroup == "Avian fecal")
Avian.group.ps.css <- prune_taxa(taxa_sums(Avian.group.ps.css) > 0, Avian.group.ps.css)
any(taxa_sums(Avian.group.ps.css)==0) # QUADRUPLE CHECKING - nope good.
Avian.group.ps.css.df <- as(sample_data(Avian.group.ps.css), "data.frame")

# Porcine 
Porcine.group.ps.css <- subset_samples(group.ps.css, SampleGroup == "Porcine fecal")
Porcine.group.ps.css <- prune_taxa(taxa_sums(Porcine.group.ps.css) > 0, Porcine.group.ps.css)
any(taxa_sums(Porcine.group.ps.css)==0) # QUADRUPLE CHECKING - nope good.
Porcine.group.ps.css.df <- as(sample_data(Porcine.group.ps.css), "data.frame")

# Human 
Human.group.ps.css <- subset_samples(group.ps.css, SampleGroup == "Human wastewater")
Human.group.ps.css <- prune_taxa(taxa_sums(Human.group.ps.css) > 0, Human.group.ps.css)
any(taxa_sums(Human.group.ps.css)==0) # QUADRUPLE CHECKING - nope good.
Human.group.ps.css.df <- as(sample_data(Human.group.ps.css), "data.frame")

# Canine 
Canine.group.ps.css <- subset_samples(group.ps.css, SampleGroup == "Canine fecal")
Canine.group.ps.css <- prune_taxa(taxa_sums(Canine.group.ps.css) > 0, Canine.group.ps.css)
any(taxa_sums(Canine.group.ps.css)==0) # QUADRUPLE CHECKING - nope good.
Canine.group.ps.css.df <- as(sample_data(Canine.group.ps.css), "data.frame")

###
####
## AllSamples - Beta diversity #######
####
### 

# Bray curtis distance matrix did converge with all samples
# Had to use a different beta diversity index, euclidean, calculated on hellinger transformed counts
AllSamples.group.ps.css.df <- as(sample_data(group.ps.css),"data.frame")
AllSamples.group.ps.css.dist <- vegdist(decostand(t(otu_table(group.ps.css)), "hell"), "euclidean") 
AllSamples.group.ps.css.ord <- vegan::metaMDS(comm = t(AllSamples.group.ps.css.dist), distance = "none", try = 100, trymax = 999, autotransform = F)

AllSamples.plt_ord_by_Group_library <- plot_ordination(group.ps.css, AllSamples.group.ps.css.ord, color = "SampleGroup", shape = "Kit") +
  theme_bw() +
  #labs(title ="Sample type") +
  stat_ellipse(aes(fill= SampleGroup), geom="polygon", alpha = 0.25) +
  geom_point(size = 2) +  # Slightly increase the point size 
  scale_color_manual(name = "Sample groups", values=c("#33a02c", "#b15928", "#6a3d9a", "#7cb9cb", "#FB9A99")) +
  scale_fill_manual(name = "Sample groups", values=c("#33a02c", "#b15928", "#6a3d9a", "#7cb9cb", "#FB9A99")) +
  theme(legend.position = "right",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        legend.text = element_text(size = 16),    # Increase legend text size
        legend.title = element_text(size = 18),   # Increase legend title size
        strip.background = element_rect(fill= "grey91", linewidth = 1.0),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
AllSamples.plt_ord_by_Group_library

# Save the final plot
ggsave("../Writing/Final_figures/Fig5_AMR_betadiversity.jpg", plot = AllSamples.plt_ord_by_Group_library, width = 20, height = 11, dpi = 300)

## Without canine
ps_sub <- subset_samples(group.ps.css, SampleGroup != "Canine fecal")
ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
PartA.group.ps.css.df <- as(sample_data(ps_sub),"data.frame")
PartA.group.ps.css.dist <- vegdist(decostand(t(otu_table(ps_sub)), "hell"), "euclidean") 
PartA.group.ps.css.ord <- vegan::metaMDS(comm = t(PartA.group.ps.css.dist), distance = "none", try = 100, trymax = 999, autotransform = F)

PartA.plt_ord_by_Group_library <- plot_ordination(group.ps.css, PartA.group.ps.css.ord, color = "SampleGroup", shape = "Kit") +
  theme_bw() +
  #labs(title ="Sample type") +
  stat_ellipse(aes(fill= SampleGroup), geom="polygon", alpha = 0.25) +
  geom_point(size = 2) +  # Slightly increase the point size 
  scale_color_manual(name = "Sample groups", values=c("#33a02c", "#b15928", "#7cb9cb", "#FB9A99")) +
  scale_fill_manual(name = "Sample groups", values=c("#33a02c", "#b15928", "#7cb9cb", "#FB9A99")) +
  theme(legend.position = "right",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        legend.text = element_text(size = 16),    # Increase legend text size
        legend.title = element_text(size = 18),   # Increase legend title size
        strip.background = element_rect(fill= "grey91", linewidth = 1.0),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
PartA.plt_ord_by_Group_library

# Save the final plot
#ggsave("../Writing/Final_figures/Fig_PartA_AMR_betadiversity.jpg", plot = PartA.plt_ord_by_Group_library, width = 20, height = 11, dpi = 300)


# permanova test of beta diversity, by SampleGroup
AllSamples.data.adonis <- adonis2(AllSamples.group.ps.css.dist ~ SampleGroup, data = group.ps.css.df, permutations = 9999)
AllSamples.data.adonis # 7e-04, R2 =0.13205 

# permanova test of beta diversity, by OriginalSample
AllSamples.data.adonis <- adonis2(AllSamples.group.ps.css.dist ~ OriginalSample, data = group.ps.css.df, permutations = 9999)
AllSamples.data.adonis # 7e-04, R2 =0.13205 



##
### Dilution ####
##
# pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
AllSamples.Dilution.permanova <- pairwise.adonis(AllSamples.group.ps.css.dist, group.ps.css.df$Dilution, perm = 9999, p.adjust.m = "BH")
AllSamples.Dilution.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
AllSamples.Dilution.disper <- betadisper(AllSamples.group.ps.css.dist, group.ps.css.df$Dilution)
AllSamples.Dilution.permdisp <- permutest(AllSamples.Dilution.disper, permutations = 9999, pairwise = T)
AllSamples.Dilution.permdisp # 

##
### InputDNA ####
##
# pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
AllSamples.InputDNA.permanova <- pairwise.adonis(AllSamples.group.ps.css.dist, group.ps.css.df$InputDNA, perm = 9999, p.adjust.m = "BH")
AllSamples.InputDNA.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
AllSamples.InputDNA.disper <- betadisper(AllSamples.group.ps.css.dist, group.ps.css.df$InputDNA)
AllSamples.InputDNA.permdisp <- permutest(AllSamples.InputDNA.disper, permutations = 9999, pairwise = T)
AllSamples.InputDNA.permdisp # looks like a few are significant

##
### Kit ####
##
# pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
AllSamples.Kit.permanova <- pairwise.adonis(AllSamples.group.ps.css.dist, group.ps.css.df$Kit, perm = 9999, p.adjust.m = "BH")
AllSamples.Kit.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
AllSamples.Kit.disper <- betadisper(AllSamples.group.ps.css.dist, group.ps.css.df$Kit)
AllSamples.Kit.permdisp <- permutest(AllSamples.Kit.disper, permutations = 9999, pairwise = T)
AllSamples.Kit.permdisp # looks like a few are significant

##
### OriginalSample ####
##
# pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
AllSamples.OriginalSample.permanova <- pairwise.adonis(AllSamples.group.ps.css.dist, group.ps.css.df$OriginalSample, perm = 9999, p.adjust.m = "BH")
AllSamples.OriginalSample.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
AllSamples.OriginalSample.disper <- betadisper(AllSamples.group.ps.css.dist, group.ps.css.df$OriginalSample)
AllSamples.OriginalSample.permdisp <- permutest(AllSamples.OriginalSample.disper, permutations = 9999, pairwise = T)
AllSamples.OriginalSample.permdisp # looks like a few are significant

##
### SampleGroup ####
##
# pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
AllSamples.SampleGroup.permanova <- pairwise.adonis(AllSamples.group.ps.css.dist, group.ps.css.df$SampleGroup, perm = 9999, p.adjust.m = "BH")
AllSamples.SampleGroup.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
AllSamples.SampleGroup.disper <- betadisper(AllSamples.group.ps.css.dist, group.ps.css.df$SampleGroup)
AllSamples.SampleGroup.permdisp <- permutest(AllSamples.SampleGroup.disper, permutations = 9999, pairwise = T)
AllSamples.SampleGroup.permdisp # looks like a few are significant

# Test with blocking
permanova_block_InputDNA <- perform_pairwise_adonis_blocking(AllSamples.group.ps.css.dist, AllSamples.group.ps.css.df, "InputDNA" ,"SampleGroup")
permanova_block_InputDNA

PERMDISP_block_InputDNA <- perform_pairwise_PERMDISP_blocking(AllSamples.group.ps.css.dist, AllSamples.group.ps.css.df, "InputDNA" ,"SampleGroup")
PERMDISP_block_InputDNA


###
####
## DatasetA - Beta diversity #######
####
### 

# Bray curtis distance matrix did converge with all samples
# Had to use a different beta diversity index, euclidean, calculated on hellinger transformed counts
DatasetA.group.ps.css.df <- as(sample_data(DatasetA.group.ps.css),"data.frame")
DatasetA.group.ps.css.dist <- vegdist(decostand(t(otu_table(DatasetA.group.ps.css)), "hell"), "euclidean") 
DatasetA.group.ps.css.ord <- vegan::metaMDS(comm = t(DatasetA.group.ps.css.dist), distance = "none", try = 100, trymax = 999, autotransform = F)

# Test with blocking - Dataset A ####
permanova_block_by_Kit <- perform_pairwise_adonis_blocking(DatasetA.group.ps.css.dist, DatasetA.group.ps.css.df, "Kit", "SampleGroup")
permanova_block_by_Kit

permdisp_block_by_Kit <- perform_pairwise_PERMDISP_blocking(DatasetA.group.ps.css.dist, DatasetA.group.ps.css.df, "Kit", "SampleGroup")
permdisp_block_by_Kit

# Test with blocking
permanova_block_by_InputDNA <- perform_pairwise_adonis_blocking(DatasetA.group.ps.css.dist, DatasetA.group.ps.css.df, "InputDNA", "SampleGroup")
permanova_block_by_InputDNA

permdisp_block_by_InputDNA <- perform_pairwise_PERMDISP_blocking(DatasetA.group.ps.css.dist, DatasetA.group.ps.css.df, "InputDNA", "SampleGroup")
permdisp_block_by_InputDNA





plt_ord_by_Dilution <- plot_ordination(DatasetA.group.ps.css, DatasetA.group.ps.css.ord, color = "SampleGroup", shape = "Kit") +
  theme_bw() +
  labs(title ="Sample type") +
  stat_ellipse(aes(fill= SampleGroup), geom="polygon", alpha = 0.25) +
  theme(legend.position = "right",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", linewidth = 1.0),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
plt_ord_by_Dilution

# permanova test of beta diversity, by SampleGroup
DatasetA.data.SampleGroup.adonis <- adonis2(DatasetA.group.ps.css.dist ~ SampleGroup, data = DatasetA.group.ps.css.df, permutations = 9999)
DatasetA.data.SampleGroup.adonis # 1e-04, R2 =0.8808


# permanova test of beta diversity, by SampleGroup
DatasetA.data.OriginalSample.adonis <- adonis2(DatasetA.group.ps.css.dist ~ OriginalSample, data = DatasetA.group.ps.css.df, permutations = 9999)
DatasetA.data.OriginalSample.adonis # 1e-04, R2 =0.9533 


# permanova test of beta diversity, by Kit
DatasetA.Kit.adonis <- adonis2(DatasetA.group.ps.css.dist ~ Kit, data = DatasetA.group.ps.css.df, permutations = 9999)
DatasetA.Kit.adonis # 


# permanova test of beta diversity, by InputDNA
DatasetA.InputDNA.adonis <- adonis2(DatasetA.group.ps.css.dist ~ InputDNA, data = DatasetA.group.ps.css.df, permutations = 9999)
DatasetA.InputDNA.adonis # 

# permanova test of beta diversity, by OriginalSample
DatasetA.OriginalSample.adonis <- adonis2(DatasetA.group.ps.css.dist ~ OriginalSample, data = DatasetA.group.ps.css.df, permutations = 9999)
DatasetA.OriginalSample.adonis #  



##
### InputDNA ####
##
# pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
DatasetA.InputDNA.permanova <- pairwise.adonis(DatasetA.group.ps.css.dist, DatasetA.group.ps.css.df$InputDNA, perm = 9999, p.adjust.m = "BH")
DatasetA.InputDNA.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
DatasetA.InputDNA.disper <- betadisper(DatasetA.group.ps.css.dist, DatasetA.group.ps.css.df$InputDNA)
DatasetA.InputDNA.permdisp <- permutest(DatasetA.InputDNA.disper, permutations = 9999, pairwise = T)
DatasetA.InputDNA.permdisp # looks like a few are significant

##
### Kit ####
##
# pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
DatasetA.Kit.permanova <- pairwise.adonis(DatasetA.group.ps.css.dist, DatasetA.group.ps.css.df$Kit, perm = 9999, p.adjust.m = "BH")
DatasetA.Kit.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
DatasetA.Kit.disper <- betadisper(DatasetA.group.ps.css.dist, DatasetA.group.ps.css.df$Kit)
DatasetA.Kit.permdisp <- permutest(DatasetA.Kit.disper, permutations = 9999, pairwise = T)
DatasetA.Kit.permdisp # looks like a few are significant

##
### OriginalSample ####
##
# pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
DatasetA.OriginalSample.permanova <- pairwise.adonis(DatasetA.group.ps.css.dist, DatasetA.group.ps.css.df$OriginalSample, perm = 9999, p.adjust.m = "BH")
DatasetA.OriginalSample.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
DatasetA.OriginalSample.disper <- betadisper(DatasetA.group.ps.css.dist, DatasetA.group.ps.css.df$OriginalSample)
DatasetA.OriginalSample.permdisp <- permutest(DatasetA.OriginalSample.disper, permutations = 9999, pairwise = T)
DatasetA.OriginalSample.permdisp # looks like a few are significant

##
### SampleGroup ####
##
# pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
DatasetA.SampleGroup.permanova <- pairwise.adonis(DatasetA.group.ps.css.dist, DatasetA.group.ps.css.df$SampleGroup, perm = 9999, p.adjust.m = "BH")
DatasetA.SampleGroup.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
DatasetA.SampleGroup.disper <- betadisper(DatasetA.group.ps.css.dist, DatasetA.group.ps.css.df$SampleGroup)
DatasetA.SampleGroup.permdisp <- permutest(DatasetA.SampleGroup.disper, permutations = 9999, pairwise = T)
DatasetA.SampleGroup.permdisp # looks like a few are significant






###
####
## Bovine - Beta diversity #######
####
### 

# Bray curtis distance matrix did converge with all samples
# Had to use a different beta diversity index, euclidean, calculated on hellinger transformed counts
Bovine.group.ps.css.dist <- vegdist(decostand(t(otu_table(Bovine.group.ps.css)), "hell"), "euclidean") 
Bovine.group.ps.css.ord <- vegan::metaMDS(comm = t(Bovine.group.ps.css.dist), distance = "none", try = 100, trymax = 999, autotransform = F)

plt_ord_by_Dilution <- plot_ordination(Bovine.group.ps.css, Bovine.group.ps.css.ord, color = "OriginalSample", shape = "Kit") +
  theme_bw() +
  labs(title ="Sample type") +
  stat_ellipse(aes(fill= OriginalSample), geom="polygon", alpha = 0.25) +
  theme(legend.position = "right",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", linewidth = 1.0),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) 
  #facet_wrap(~ InputDNA)
plt_ord_by_Dilution

# permanova test of beta diversity, by Kit
Bovine.Kit.adonis <- adonis2(Bovine.group.ps.css.dist ~ Kit, data = Bovine.group.ps.css.df, permutations = 9999)
  Bovine.Kit.adonis # 

# permanova test of beta diversity, by Dilution
Bovine.Dilution.adonis <- adonis2(Bovine.group.ps.css.dist ~ Dilution, data = Bovine.group.ps.css.df, permutations = 9999)
Bovine.Dilution.adonis #  

# permanova test of beta diversity, by InputDNA
Bovine.InputDNA.adonis <- adonis2(Bovine.group.ps.css.dist ~ InputDNA, data = Bovine.group.ps.css.df, permutations = 9999)
Bovine.InputDNA.adonis # 

# permanova test of beta diversity, by OriginalSample
Bovine.OriginalSample.adonis <- adonis2(Bovine.group.ps.css.dist ~ OriginalSample, data = Bovine.group.ps.css.df, permutations = 9999)
Bovine.OriginalSample.adonis #  


##
### InputDNA ####
##
# pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
InputDNA.Bovine.permanova <- pairwise.adonis(Bovine.group.ps.css.dist, Bovine.group.ps.css.df$InputDNA, perm = 9999, p.adjust.m = "BH")
InputDNA.Bovine.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
InputDNA.Bovine.disper <- betadisper(Bovine.group.ps.css.dist, Bovine.group.ps.css.df$InputDNA)
InputDNA.Bovine.permdisp <- permutest(InputDNA.Bovine.disper, permutations = 9999, pairwise = T)
InputDNA.Bovine.permdisp # looks like a few are significant

##
### Kit ####
##
# pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
Kit.Bovine.permanova <- pairwise.adonis(Bovine.group.ps.css.dist, Bovine.group.ps.css.df$Kit, perm = 9999, p.adjust.m = "BH")
Kit.Bovine.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
Kit.Bovine.disper <- betadisper(Bovine.group.ps.css.dist, Bovine.group.ps.css.df$Kit)
Kit.Bovine.permdisp <- permutest(Kit.Bovine.disper, permutations = 9999, pairwise = T)
Kit.Bovine.permdisp # looks like a few are significant

##
### OriginalSample ####
##
# pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
OriginalSample.Bovine.permanova <- pairwise.adonis(Bovine.group.ps.css.dist, Bovine.group.ps.css.df$OriginalSample, perm = 9999, p.adjust.m = "BH")
OriginalSample.Bovine.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
OriginalSample.Bovine.disper <- betadisper(Bovine.group.ps.css.dist, Bovine.group.ps.css.df$OriginalSample)
OriginalSample.Bovine.permdisp <- permutest(OriginalSample.Bovine.disper, permutations = 9999, pairwise = T)
OriginalSample.Bovine.permdisp # looks like a few are significant



###
####
## Avian - Beta diversity #######
####
###


# Bray curtis distance matrix did converge with all samples
# Had to use a different beta diversity index, euclidean, calculated on hellinger transformed counts
Avian.group.ps.css.dist <- vegdist(decostand(t(otu_table(Avian.group.ps.css)), "hell"), "euclidean") 
Avian.group.ps.css.ord <- vegan::metaMDS(comm = t(Avian.group.ps.css.dist), distance = "none", try = 100, trymax = 999, autotransform = F)

plt_ord_by_Dilution <- plot_ordination(Avian.group.ps.css, Avian.group.ps.css.ord, color = "OriginalSample", shape = "Kit") +
  theme_bw() +
  labs(title ="Sample type") +
  stat_ellipse(aes(fill= OriginalSample), geom="polygon", alpha = 0.25) +
  theme(legend.position = "right",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", linewidth = 1.0),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) 
#facet_wrap(~ InputDNA)
plt_ord_by_Dilution


# permanova test of beta diversity, by Kit
Avian.Kit.adonis <- adonis2(Avian.group.ps.css.dist ~ Kit, data = Avian.group.ps.css.df, permutations = 9999)
Avian.Kit.adonis # 


# permanova test of beta diversity, by InputDNA
Avian.InputDNA.adonis <- adonis2(Avian.group.ps.css.dist ~ InputDNA, data = Avian.group.ps.css.df, permutations = 9999)
Avian.InputDNA.adonis # 

# permanova test of beta diversity, by OriginalSample
Avian.OriginalSample.adonis <- adonis2(Avian.group.ps.css.dist ~ OriginalSample, data = Avian.group.ps.css.df, permutations = 9999)
Avian.OriginalSample.adonis #


##
### InputDNA ####
##
# pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
InputDNA.Avian.permanova <- pairwise.adonis(Avian.group.ps.css.dist, Avian.group.ps.css.df$InputDNA, perm = 9999, p.adjust.m = "BH")
InputDNA.Avian.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
InputDNA.Avian.disper <- betadisper(Avian.group.ps.css.dist, Avian.group.ps.css.df$InputDNA)
InputDNA.Avian.permdisp <- permutest(InputDNA.Avian.disper, permutations = 9999, pairwise = T)
InputDNA.Avian.permdisp # looks like a few are significant

##
### Kit ####
##
# pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
Kit.Avian.permanova <- pairwise.adonis(Avian.group.ps.css.dist, Avian.group.ps.css.df$Kit, perm = 9999, p.adjust.m = "BH")
Kit.Avian.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
Kit.Avian.disper <- betadisper(Avian.group.ps.css.dist, Avian.group.ps.css.df$Kit)
Kit.Avian.permdisp <- permutest(Kit.Avian.disper, permutations = 9999, pairwise = T)
Kit.Avian.permdisp # looks like a few are significant

##
### OriginalSample ####
##
# pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
OriginalSample.Avian.permanova <- pairwise.adonis(Avian.group.ps.css.dist, Avian.group.ps.css.df$OriginalSample, perm = 9999, p.adjust.m = "BH")
OriginalSample.Avian.permanova # there are significant differences between some groups


## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
OriginalSample.Avian.disper <- betadisper(Avian.group.ps.css.dist, Avian.group.ps.css.df$OriginalSample)
OriginalSample.Avian.permdisp <- permutest(OriginalSample.Avian.disper, permutations = 9999, pairwise = T)
OriginalSample.Avian.permdisp # looks like a few are significant

###
####
## Porcine - Beta diversity #######
####
###


# Bray curtis distance matrix did converge with all samples
# Had to use a different beta diversity index, euclidean, calculated on hellinger transformed counts
Porcine.group.ps.css.dist <- vegdist(decostand(t(otu_table(Porcine.group.ps.css)), "hell"), "euclidean") 
Porcine.group.ps.css.ord <- vegan::metaMDS(comm = t(Porcine.group.ps.css.dist), distance = "none", try = 100, trymax = 999, autotransform = F)

plt_ord_by_Dilution <- plot_ordination(Porcine.group.ps.css, Porcine.group.ps.css.ord, color = "OriginalSample", shape = "Kit") +
  theme_bw() +
  labs(title ="Sample type") +
  stat_ellipse(aes(fill= OriginalSample), geom="polygon", alpha = 0.25) +
  theme(legend.position = "right",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", linewidth = 1.0),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) 
#facet_wrap(~ InputDNA)
plt_ord_by_Dilution

# permanova test of beta diversity, by Dilution
data.adonis <- adonis2(Porcine.group.ps.css.dist ~ Dilution, data = Porcine.group.ps.css.df, permutations = 9999)
data.adonis # 7e-04, R2 =0.13205 


# permanova test of beta diversity, by Kit
Porcine.Kit.adonis <- adonis2(Porcine.group.ps.css.dist ~ Kit, data = Porcine.group.ps.css.df, permutations = 9999)
Porcine.Kit.adonis # 


# permanova test of beta diversity, by InputDNA
Porcine.InputDNA.adonis <- adonis2(Porcine.group.ps.css.dist ~ InputDNA, data = Porcine.group.ps.css.df, permutations = 9999)
Porcine.InputDNA.adonis # 

# permanova test of beta diversity, by OriginalSample
Porcine.OriginalSample.adonis <- adonis2(Porcine.group.ps.css.dist ~ OriginalSample, data = Porcine.group.ps.css.df, permutations = 9999)
Porcine.OriginalSample.adonis #

##
### Dilution ####
##
# pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
Dilution.Porcine.permanova <- pairwise.adonis(Porcine.group.ps.css.dist, Porcine.group.ps.css.df$Dilution, perm = 9999, p.adjust.m = "BH")
Dilution.Porcine.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
Dilution.Porcine.disper <- betadisper(Porcine.group.ps.css.dist, Porcine.group.ps.css.df$Dilution)
Dilution.Porcine.permdisp <- permutest(Dilution.Porcine.disper, permutations = 9999, pairwise = T)
Dilution.Porcine.permdisp # looks like a few are significant

##
### InputDNA ####
##
# pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
InputDNA.Porcine.permanova <- pairwise.adonis(Porcine.group.ps.css.dist, Porcine.group.ps.css.df$InputDNA, perm = 9999, p.adjust.m = "BH")
InputDNA.Porcine.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
InputDNA.Porcine.disper <- betadisper(Porcine.group.ps.css.dist, Porcine.group.ps.css.df$InputDNA)
InputDNA.Porcine.permdisp <- permutest(InputDNA.Porcine.disper, permutations = 9999, pairwise = T)
InputDNA.Porcine.permdisp # looks like a few are significant

##
### Kit ####
##
# pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
Kit.Porcine.permanova <- pairwise.adonis(Porcine.group.ps.css.dist, Porcine.group.ps.css.df$Kit, perm = 9999, p.adjust.m = "BH")
Kit.Porcine.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
Kit.Porcine.disper <- betadisper(Porcine.group.ps.css.dist, Porcine.group.ps.css.df$Kit)
Kit.Porcine.permdisp <- permutest(Kit.Porcine.disper, permutations = 9999, pairwise = T)
Kit.Porcine.permdisp # looks like a few are significant

##
### OriginalSample ####
##
# pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
OriginalSample.Porcine.permanova <- pairwise.adonis(Porcine.group.ps.css.dist, Porcine.group.ps.css.df$OriginalSample, perm = 9999, p.adjust.m = "BH")
OriginalSample.Porcine.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
OriginalSample.Porcine.disper <- betadisper(Porcine.group.ps.css.dist, Porcine.group.ps.css.df$OriginalSample)
OriginalSample.Porcine.permdisp <- permutest(OriginalSample.Porcine.disper, permutations = 9999, pairwise = T)
OriginalSample.Porcine.permdisp # looks like a few are significant

###
####
## Human - Beta diversity #######
####
###


# Bray curtis distance matrix did converge with all samples
# Had to use a different beta diversity index, euclidean, calculated on hellinger transformed counts
Human.group.ps.css.dist <- vegdist(decostand(t(otu_table(Human.group.ps.css)), "hell"), "euclidean") 
Human.group.ps.css.ord <- vegan::metaMDS(comm = t(Human.group.ps.css.dist), distance = "none", try = 100, trymax = 999, autotransform = F)

plt_ord_by_Dilution <- plot_ordination(Human.group.ps.css, Human.group.ps.css.ord, color = "OriginalSample", shape = "Kit") +
  theme_bw() +
  labs(title ="Sample type") +
  stat_ellipse(aes(fill= OriginalSample), geom="polygon", alpha = 0.25) +
  theme(legend.position = "right",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", linewidth = 1.0),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) 
#facet_wrap(~ InputDNA)
plt_ord_by_Dilution

# permanova test of beta diversity, by Dilution
data.adonis <- adonis2(Human.group.ps.css.dist ~ Dilution, data = Human.group.ps.css.df, permutations = 9999)
data.adonis # 7e-04, R2 =0.13205 

# permanova test of beta diversity, by Kit
Human.Kit.adonis <- adonis2(Human.group.ps.css.dist ~ Kit, data = Human.group.ps.css.df, permutations = 9999)
Human.Kit.adonis # 


# permanova test of beta diversity, by InputDNA
Human.InputDNA.adonis <- adonis2(Human.group.ps.css.dist ~ InputDNA, data = Human.group.ps.css.df, permutations = 9999)
Human.InputDNA.adonis # 

# permanova test of beta diversity, by OriginalSample
Human.OriginalSample.adonis <- adonis2(Human.group.ps.css.dist ~ OriginalSample, data = Human.group.ps.css.df, permutations = 9999)
Human.OriginalSample.adonis #

##
### InputDNA ####
##
# pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
InputDNA.Human.permanova <- pairwise.adonis(Human.group.ps.css.dist, Human.group.ps.css.df$InputDNA, perm = 9999, p.adjust.m = "BH")
InputDNA.Human.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
InputDNA.Human.disper <- betadisper(Human.group.ps.css.dist, Human.group.ps.css.df$InputDNA)
InputDNA.Human.permdisp <- permutest(InputDNA.Human.disper, permutations = 9999, pairwise = T)
InputDNA.Human.permdisp # looks like a few are significant

##
### Kit ####
##
# pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
Kit.Human.permanova <- pairwise.adonis(Human.group.ps.css.dist, Human.group.ps.css.df$Kit, perm = 9999, p.adjust.m = "BH")
Kit.Human.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
Kit.Human.disper <- betadisper(Human.group.ps.css.dist, Human.group.ps.css.df$Kit)
Kit.Human.permdisp <- permutest(Kit.Human.disper, permutations = 9999, pairwise = T)
Kit.Human.permdisp # looks like a few are significant

##
### OriginalSample ####
##
# pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
OriginalSample.Human.permanova <- pairwise.adonis(Human.group.ps.css.dist, Human.group.ps.css.df$OriginalSample, perm = 9999, p.adjust.m = "BH")
OriginalSample.Human.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
OriginalSample.Human.disper <- betadisper(Human.group.ps.css.dist, Human.group.ps.css.df$OriginalSample)
OriginalSample.Human.permdisp <- permutest(OriginalSample.Human.disper, permutations = 9999, pairwise = T)
OriginalSample.Human.permdisp # looks like a few are significant

###
####
## Canine - Beta diversity #######
####
###


# Bray curtis distance matrix did converge with all samples
# Had to use a different beta diversity index, euclidean, calculated on hellinger transformed counts
Canine.group.ps.css.dist <- vegdist(decostand(t(otu_table(Canine.group.ps.css)), "hell"), "euclidean") 
Canine.group.ps.css.ord <- vegan::metaMDS(comm = t(Canine.group.ps.css.dist), distance = "none", try = 100, trymax = 999, autotransform = F)

plt_ord_by_Dilution <- plot_ordination(Canine.group.ps.css, Canine.group.ps.css.ord, color = "OriginalSample", shape = "Kit") +
  theme_bw() +
  labs(title ="Sample type") +
  stat_ellipse(aes(fill= OriginalSample), geom="polygon", alpha = 0.25) +
  theme(legend.position = "right",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", linewidth = 1.0),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) 
#facet_wrap(~ InputDNA)
plt_ord_by_Dilution

# permanova test of beta diversity, by Dilution
data.adonis <- adonis2(Canine.group.ps.css.dist ~ InputDNA, data = Canine.group.ps.css.df, permutations = 9999)
data.adonis # 7e-04, R2 =0.13205 

# permanova test of beta diversity, by Dilution
data.adonis <- adonis2(Canine.group.ps.css.dist ~ OriginalSample, data = Canine.group.ps.css.df, permutations = 9999)
data.adonis # 

##
### InputDNA ####
##
# pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
InputDNA.Canine.permanova <- pairwise.adonis(Canine.group.ps.css.dist, Canine.group.ps.css.df$InputDNA, perm = 9999, p.adjust.m = "BH")
InputDNA.Canine.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
InputDNA.Canine.disper <- betadisper(Canine.group.ps.css.dist, Canine.group.ps.css.df$InputDNA)
InputDNA.Canine.permdisp <- permutest(InputDNA.Canine.disper, permutations = 9999, pairwise = T)
InputDNA.Canine.permdisp # looks like a few are significant
