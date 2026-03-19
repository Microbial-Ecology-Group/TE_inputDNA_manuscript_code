library(ggsci) # if not using this, make sure to switch out the palatte in the figures

#############################################################################################
####################################   CSS TRANSFORM    ###################################
sample_data(data_noSNP)$Dilution <- factor(sample_data(data_noSNP)$Dilution, levels = c("Standard","1:2 dilution",
                                                                                        "1:3 dilution", "1:4 dilution","1:5 dilution","1:10 dilution"))

data_noSNP <- prune_taxa(taxa_sums(data_noSNP) > 0, data_noSNP)
any(taxa_sums(data_noSNP)==0) # QUADRUPLE CHECKING - nope good.

data_noSNP.css <- phyloseq_transform_css(data_noSNP, log = F)
data_noSNP.css.df <- as(sample_data(data_noSNP.css), "data.frame")



#############################################################################################
##############################         RELA ABUNDANCE         ###############################
#############################################################################################
#############################################################################################

rel_abund <- transform_sample_counts(data_noSNP.css, function(x) {x/sum(x)}*100)

  # agglomerate
ra_class <- tax_glom(rel_abund, taxrank = "class")
ra_class # 40 classes, 153 samples


ra_class_melt_temp <- psmelt(ra_class)


ps_class.css <- tax_glom(data_noSNP.css, taxrank = "class")
ps_class.melt <- psmelt(ps_class.css)
#ra_class_melt_temp <- ps_class.melt

# ra_mechanism <- tax_glom(rel_abund, taxrank = "Mechanism")
# ra_mechanism # 131 mechanisms
# ra_mechanism_melt <- psmelt(ra_mechanism)
# 
# ra_group <- tax_glom(rel_abund, taxrank = "Group")
# ra_group # 607 groups
# ra_group_melt <- psmelt(ra_group)



######## Re-label low abundance taxa into single category

# convert Phylum to a character vector from a factor because R
ra_class_melt_temp$class <- as.character(ra_class_melt_temp$class)

# calculate percentage across all counts
all_sample_factor_by_abund <- ra_class_melt_temp %>%
  dplyr::group_by(class) %>%
  dplyr::summarize(median_class = median(Abundance)) %>%
  arrange(-median_class)

# find Phyla whose rel. abund. is less than 1%
remainder <- all_sample_factor_by_abund[all_sample_factor_by_abund$median_class <= 0.3,]$class

# change the name of taxa to whose rel. abund. is less than 1% to "Low abundance phyla (<1%)"
ra_class_melt_temp[ra_class_melt_temp$class %in% remainder,]$class <- 'Low abundance class (<3%)'

# Determine the number of taxa
length(unique(ra_class_melt_temp$class))

# ra_class_melt <- ps_class.melt %>%
#   group_by(Sample, Class) %>%
#   mutate(Rel_Abundance = Abundance/  sum(Abundance))



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

# This requires a directory called figures to be in your working directory
# All 3 lines have to be run for the file to be created.
png("Figures/Allsamples_NMDS_Dilution.png", width = 1200, height = 800)
plt_ord_by_Dilution
dev.off()

# svg format is nice because it stores all graphic data and can be requested by some manuscripts
svg("Figures/Allsamples_NMDS_Dilution.svg", width = 500, height = 500)
plt_ord_by_Dilution
dev.off()



##
### Clustering at ASV level ####
##

# Cluster using the function hclust() and default settings
ps_AMR.dist <- vegdist(t(otu_table(data_noSNP.css)), method = "bray")
ps_AMR.hclust <- hclust(ps_AMR.dist)
plot(ps_AMR.hclust) # example plot

# Extract data as dendrogram
ps_AMR.dendro <- as.dendrogram(ps_AMR.hclust)
ps_AMR.dendro.data <- dendro_data(ps_AMR.dendro, type = "rectangle")
ps_AMR.dendro.data #  this object contains $segments and $labels

# Sample names in order based on clustering
sample_names_ps_AMR <- ps_AMR.dendro.data$labels$label

# Add metadata
ps_AMR.dendro.data$labels <- left_join(ps_AMR.dendro.data$labels, sample_data(data_noSNP.css.df), by = c("label" = "Sample"))

# Setup the data, so that the layout is inverted (this is more 
# "clear" than simply using coord_flip())
segment_data_ps_AMR <- with(
  segment(ps_AMR.dendro.data),
  data.frame(x = y, y = x, xend = yend, yend = xend))

# Use the dendrogram label data to position the gene labels
gene_pos_table_ps_AMR <- with(
  ps_AMR.dendro.data$labels, 
  data.frame(y_center = x, gene = as.character(label), x = y , trt = as.character(Dilution), treatment_order = as.character(Dilution) , height = 1))

# Table to position the samples
sample_pos_table_ps_AMR <- data.frame(sample = sample_names_ps_AMR) %>%
  dplyr::mutate(x_center = (1:dplyr::n()),  width = 1)

##
######## Relative abundance bar plot #########
##

ra_class_melt <- ra_class_melt_temp

# ra_class_melt <- ra_class_melt_temp %>%
#   group_by(Class, Sample) %>%
#   mutate(Abundance = sum(Abundance))

# Use class melted data and add gene locations
joined_ra_class_ps_class_melt <- ra_class_melt %>%
  left_join(gene_pos_table_ps_AMR, by = c("Sample" = "gene")) %>%
  left_join(sample_pos_table_ps_AMR, by = c("Sample" = "sample")) 

# Calculate the mean relative abundance of each class taxa, sort by most abundant to least

ra_class_melt$class <- as.character(ra_class_melt$class)

factor_by_abund <- ra_class_melt %>%
  dplyr::group_by(class) %>%
  dplyr::summarize(median_class = median(Abundance)) %>%
  arrange(-median_class)

# Use the sorted class levels to clean up the order of taxa in the relative abundance plots
joined_ra_class_ps_class_melt$class <- factor(joined_ra_class_ps_class_melt$class, levels = as.character(factor_by_abund$class))

# Limits for the vertical axes
gene_axis_limits_ps_class <- with(
  gene_pos_table_ps_AMR, 
  c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))) +   0.1 * c(-1, 1) # extra spacing: 0.1


# Useful 20 colors + black and white (https://sashamaps.net/docs/resources/20-colors/)
#'#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#ffffff', '#000000'

#
## Relative abundance plot
#
plt_rel_class_ps_class <- ggplot(joined_ra_class_ps_class_melt, 
                                          aes(x = x_center, y = Abundance, fill = class, 
                                              height = height, width = width)) + 
  coord_flip() +
  geom_bar(stat = "identity", colour = "black") +
  #scale_fill_manual(values = col_vector) + #use this if not using color palette from ggthemes of ggsci
  scale_fill_tableau(palette = "Tableau 20")+
  scale_x_continuous(breaks = sample_pos_table_ps_AMR$x_center, 
                     labels = sample_pos_table_ps_AMR$sample, 
                     expand = c(0, 0)) + 
  # For the y axis, alternatively set the labels as: gene_position_table$gene
  #scale_y_continuous(breaks = gene_pos_table[, "y_center"], labels = rep("", nrow(gene_pos_table)),limits = gene_axis_limits, expand = c(0, 0)) + 
  labs(x = "", y = "Relative abundance") +
  theme_bw() +
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.75),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.75, lineend = "square", colour = "black"),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, colour = "black"),
        plot.margin = unit(c(1, 0.2, 0.2, -0.7), "cm"), 
        panel.grid.minor = element_blank())
plt_rel_class_ps_class

# Dendrogram plot
plt_dendr_class_ps_class <- ggplot(segment_data_ps_AMR) + 
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend), size =.75, lineend = "round", linejoin = "round") +
  geom_point(data = gene_pos_table_ps_AMR, aes(x,y_center, colour = trt, fill = trt, shape = treatment_order),
             size = 8, stroke = 0.75, position = position_nudge(x = -0.03)) +
  scale_y_continuous(breaks = gene_pos_table_ps_AMR$y_center, 
                     labels = gene_pos_table_ps_AMR$gene, 
                     limits = gene_axis_limits_ps_class, 
                     expand = c(0, 0)) + 
  labs(x = "Ward's Distance", y = "", colour = "", size = "", title = "Sample type") +
  scale_x_reverse() + 
  scale_color_manual(values=c("orange","red",  "aquamarine4","blue","brown","yellow")) +
  scale_shape_manual(values =c(15,15,15,15,15,15)) +
  theme_bw() + 
  theme(legend.position = "right",
        panel.border = element_blank(),
        plot.title = element_text(size = 30),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(size = 0.75),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, colour = "black"))
plt_dendr_class_ps_class

class_ps_class_plots <- plot_grid(plt_dendr_class_ps_class, plt_rel_class_ps_class, align = 'h', rel_widths = c(0.5, 1.5))

# This requires a directory called figures to be in your working directory
# All 3 lines have to be run for the file to be created.
png("Figures/AMR_Class_AbunDendro_byDilution_fig.png", width = 1200, height = 800)
class_ps_class_plots
dev.off()

# svg format is nice because it stores all graphic data and can be requested by some manuscripts
svg("Figures/Feces_Class_EndWashout_byDilution_fig.svg", width = 500, height = 500)
class_ps_class_plots
dev.off()


