# Updated relative abundance plot template

library(ggsci) # if not using this, make sure to switch out the palatte in the figures
library(metagMisc)
library(microViz)
library(ggdendro)

all_AMR_group.ps <- merge_phyloseq(shotgun_AMR_group.ps,group.ps)
no_canine_group.ps <- subset_samples(
  all_AMR_group.ps,
  SampleGroup != "Canine fecal" )

no_canine_group.ps <- prune_taxa(taxa_sums(no_canine_group.ps) > 0, no_canine_group.ps)



sample_data(no_canine_group.ps)$SampleGroup

my_palette <- c( "#B2DF8A", "#33A02C","#A6CEE3", "#1F78B4", "#E31A1C", "#FDBF6F", 
                 "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#1ff8ff", "#1B9E77", 
                 "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", 
                 "#4b6a53", "#b249d5", "#7edc45", "#5c47b8", "#cfd251", "#ff69b4", "#69c86c", 
                 "#cd3e50", "#83d5af", "#da6130", "#5e79b2", "#c29545", "#532a5a", "#5f7b35", "#FB9A99",
                 "#c497cf", "#773a27", "#7cb9cb", "#594e50", "#d3c4a8", "#c17e7f", "lightgrey")



# palettes -----------------------------------------------------
dna_pal  <- c(`Shotgun` = "#00007b", `200ng` = "#ffd700",
              `400ng`   = "#ffb74c", `800ng` = "#e26700")

sg_pal   <- c(`Avian fecal`        = "#33a02c",
              `Bovine fecal`       = "#8c510a",
              #`Canine fecal`       = "#6a3d9a",
              `Human wastewater`   = "#67a9cf",
              `Porcine fecal`      = "#ef8a62")

type_pal <- c(`Shotgun` = "#00007b", `TE` = "#ffd700")  # same hues as dna


#############################################################################################
####################################   CSS TRANSFORM    ###################################



# Check the minimum sample sum first
mean(sample_sums(no_canine_group.ps))
min(sample_sums(no_canine_group.ps))
min_depth <- min(sample_sums(no_canine_group.ps))
cat("Rarefying to:", min_depth, "reads\n")

# Rarefy to minimum depth
no_canine_group.ps.rare <- rarefy_even_depth(
  no_canine_group.ps,
  sample.size = min_depth,
  rngseed = 42,        # for reproducibility
  replace = FALSE,     # sampling without replacement
  trimOTUs = TRUE,     # remove OTUs that disappear after rarefaction
  verbose = TRUE
)

# Extract sample data as data frame (same as your CSS workflow)
no_canine_group.ps.rare.df <- as(sample_data(no_canine_group.ps.rare), "data.frame")


#no_canine_group.ps.css <- phyloseq_transform_css(no_canine_group.ps, log = F)
#no_canine_group.ps.css.df <- as(sample_data(no_canine_group.ps.css), "data.frame")


#############################################################################################
##############################         RELA ABUNDANCE         ###############################
#############################################################################################
#############################################################################################

rel_abund <- transform_sample_counts(no_canine_group.ps.rare, function(x) {x/sum(x)}*100)

# agglomerate
ra_class <- tax_glom(rel_abund, taxrank = "class")
ra_class # 52 taxa and 110 samples, rarefied 48 taxa


ra_class_melt_temp <- psmelt(ra_class)


ps_class.css <- tax_glom(no_canine_group.ps.rare, taxrank = "class")
ps_class.css <- prune_taxa(taxa_sums(ps_class.css) > 1, ps_class.css)


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

# Convert 'class' to a character vector from a factor to avoid issues with factor levels
ra_class_melt_temp$class <- as.character(ra_class_melt_temp$class)

# Step 1: Calculate median abundance for each class across all samples
all_sample_factor_by_abund <- ra_class_melt_temp %>%
  dplyr::group_by(class) %>%
  dplyr::summarize(median_class = median(Abundance)) %>%
  arrange(-median_class)

# Step 2: Identify classes where the median relative abundance is <= 1%
remainder <- all_sample_factor_by_abund %>%
  dplyr::filter(median_class <= 1) %>%
  dplyr::pull(class)

# Step 3: Relabel low-abundance classes as 'Low abundance classes (<1%)'
ra_class_melt_temp$class <- ifelse(ra_class_melt_temp$class %in% remainder, 
                                   'Low abundance classes (<1%)', 
                                   ra_class_melt_temp$class)

# Step 4: Group by Sample and relabeled 'class', sum Abundance (already relative abundance, so no recalculation)
ra_class_melt_temp <- ra_class_melt_temp %>%
  dplyr::group_by(Sample, class) %>%
  dplyr::summarize(Abundance = sum(Abundance), .groups = "drop")

# View the resulting data frame
length(unique(ra_class_melt_temp$class))


### test


##
### Clustering at ASV level ####
##

# Cluster using the function hclust() and default settings
ps_AMR.dist <- vegdist(t(otu_table(no_canine_group.ps.rare)), method = "bray")
ps_AMR.hclust <- hclust(ps_AMR.dist)
plot(ps_AMR.hclust) # example plot

# Extract data as dendrogram
ps_AMR.dendro <- as.dendrogram(ps_AMR.hclust)
ps_AMR.dendro.data <- dendro_data(ps_AMR.dendro, type = "rectangle")
ps_AMR.dendro.data #  this object contains $segments and $labels

# Sample names in order based on clustering
sample_names_ps_AMR <- ps_AMR.dendro.data$labels$label

# Add metadata
ps_AMR.dendro.data$labels <- left_join(ps_AMR.dendro.data$labels, sample_data(no_canine_group.ps.rare.df), by = c("label" = "Sample"))

# Setup the data, so that the layout is inverted (this is more 
# "clear" than simply using coord_flip())
segment_data_ps_AMR <- with(
  segment(ps_AMR.dendro.data),
  data.frame(x = y, y = x, xend = yend, yend = xend))

# Use the dendrogram label data to position the gene labels
gene_pos_table_ps_AMR <- with(
  ps_AMR.dendro.data$labels, 
  data.frame(y_center = x, gene = as.character(label), x = y , trt = as.character(SampleGroup),
             Type = as.character(Type) , 
             InputDNA = as.character(InputDNA),
             height = 1))

gene_pos_table_ps_AMR
gene_pos_table_ps_AMR$trt <- factor(gene_pos_table_ps_AMR$trt, levels = c("Human wastewater","Porcine fecal", "Avian fecal", "Bovine fecal","Canine fecal"))
gene_pos_table_ps_AMR$Type <- factor(gene_pos_table_ps_AMR$Type, levels = c("Shotgun","TE"))
gene_pos_table_ps_AMR$InputDNA <- factor(gene_pos_table_ps_AMR$InputDNA, levels = c("Shotgun","200ng","400ng","800ng"))


# Table to position the samples
sample_pos_table_ps_AMR <- data.frame(sample = sample_names_ps_AMR) %>%
  dplyr::mutate(x_center = (1:dplyr::n()),  width = 1)

##
######## Relative abundance bar plot #########
##

ra_class_melt <- ra_class_melt_temp

# ra_class_melt <- ra_class_melt_temp %>%
#   dplyr::group_by(class, Sample) %>%
#   dplyr::mutate(Abundance = sum(Abundance))

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


# Reorder class factors

# Step 1: Calculate the median Abundance for each class across all samples (excluding 'Low abundance classes (<1%)')
median_abundance_by_class <- joined_ra_class_ps_class_melt %>%
  dplyr::filter(class != "Low abundance classes (<1%)") %>%
  dplyr::group_by(class) %>%
  dplyr::summarize(median_abundance = median(Abundance)) %>%
  arrange(-median_abundance)  # Sort in descending order of median abundance

# Step 2: Check the actual class names in 'median_abundance_by_class'
print(median_abundance_by_class$class)

# Step 3: Capture the class names in the desired order
# Add 'Low abundance classes (<1%)' at the end of the levels
all_levels <- c(as.character(median_abundance_by_class$class), "Low abundance classes (<1%)")

# Step 4: Ensure that the 'class' column has the correct levels
# Reorder the factor levels of 'class' based on the computed order
joined_ra_class_ps_class_melt$class <- factor(joined_ra_class_ps_class_melt$class,
                                              levels = all_levels,   # Use class names explicitly
                                              ordered = TRUE)

# Step 5: Verify if all levels are retained correctly
print(levels(joined_ra_class_ps_class_melt$class))  # This should show all levels including 'Low abundance classes (<1%)'

class_names <- as.character(levels(joined_ra_class_ps_class_melt$class))  # This should show all levels including 'Low abundance classes (<1%)'



named_palette <- setNames(my_palette[1:length(class_names)], class_names)

# Print the class names and their corresponding colors
named_palette
named_palette["Low abundance classes (<1%)"] <- "gray50"



#
## Relative abundance plot
#
plt_rel_class_ps_class <- ggplot(joined_ra_class_ps_class_melt, 
                                 aes(x = x_center, y = Abundance, fill = class, 
                                     height = height, width = width)) + 
  coord_flip() +
  geom_bar(stat = "identity", colour = "black", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = named_palette, name = "AMR class") + #use this if not using color palette from ggthemes of ggsci
  #scale_fill_tableau(palette = "Tableau 20",name = "AMR class")+
  scale_x_continuous(breaks = sample_pos_table_ps_AMR$x_center, 
                     labels = sample_pos_table_ps_AMR$sample, 
                     expand = c(0, 0)) + 
  labs(x = "", y = "Relative abundance") +
  theme_bw() +
  theme(legend.position = "right",
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22),
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


## Dendrogram plot ####

### second option ####
library(ggnewscale)      # ➊ lets us add “new” fill / colour scales
# --------------------------------------------------------------

# --------------------------------------------------------------
bar_w <- 0.04
gene_pos_table_ps_AMR <- gene_pos_table_ps_AMR %>%
  dplyr::mutate(x_dna   = x - 0.05,
                x_group = x + 0.05,
                x_type  = x )


plt_dendr_class_ps_class <- ggplot(segment_data_ps_AMR) +
  ## dendrogram branches
  geom_segment(aes(x, y, xend = xend, yend = yend),
               linewidth = .75, lineend = "round", linejoin = "round")

# ───────────────────────── (1) Input-DNA tiles ─────────────────────────
plt_dendr_class_ps_class <- plt_dendr_class_ps_class +
  geom_tile(data   = gene_pos_table_ps_AMR,
            aes(x = x_dna, y = y_center, fill = InputDNA),
            width  = bar_w, height = 1, colour = NA) +
  scale_fill_manual("Input DNA (ng)", values = dna_pal, drop = FALSE,guide  = guide_legend(order = 3)) +
  new_scale_fill()                    # ➋ start a *new* fill scale

# ───────────────────────── (2) Sample-group tiles ──────────────────────
plt_dendr_class_ps_class <- plt_dendr_class_ps_class +
  geom_tile(data   = gene_pos_table_ps_AMR,
            aes(x = x_group, y = y_center, fill = trt),
            width  = bar_w, height = 1, colour = NA) +
  scale_fill_manual("Sample group", values = sg_pal, drop = FALSE, guide  = guide_legend(order = 1)) +
  new_scale_fill()

# ───────────────────────── (3) Sequencing-type tiles ───────────────────
plt_dendr_class_ps_class <- plt_dendr_class_ps_class +
  geom_tile(data   = gene_pos_table_ps_AMR,
            aes(x = x_type, y = y_center, fill = Type),
            width  = bar_w, height = 1, colour = NA) +
  scale_fill_manual("Sequencing type", values = type_pal, drop = FALSE,guide  = guide_legend(order = 2))

# ───────────────────────── axes + theme ────────────────────────────────
plt_dendr_class_ps_class <- plt_dendr_class_ps_class +
  scale_y_continuous(breaks = gene_pos_table_ps_AMR$y_center,
                     labels = gene_pos_table_ps_AMR$gene,
                     limits = gene_axis_limits_ps_class,
                     expand = c(0, 0)) +
  scale_x_reverse() +
  labs(x = "Bray-Curtis distance", y = NULL) +
  theme_bw() +
  theme(panel.border   = element_blank(),
        panel.grid     = element_blank(),
        axis.text.y    = element_blank(),
        axis.ticks.y   = element_blank(),
        axis.line.x    = element_line(linewidth=.75),
        legend.title   = element_text(size = 22),
        legend.text    = element_text(size = 18),
        legend.position = "right")

print(plt_dendr_class_ps_class)


## Combine relative abundance plot and dendrogram

# Extract legends
legend1 <- cowplot::get_legend(plt_dendr_class_ps_class)
legend2 <- cowplot::get_legend(plt_rel_class_ps_class)

# Remove legends from the original plots
plt_dendr_class_ps_class <- plt_dendr_class_ps_class + theme(legend.position = "none")
plt_rel_class_ps_class <- plt_rel_class_ps_class + theme(legend.position = "none")

## make empty plots to help with spacing legends
empty_plot1 <- ggplot() + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
empty_plot2 <- ggplot() + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
empty_plot3 <- ggplot() + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))

# Arrange legends above one another with spacers
legends <- plot_grid(empty_plot1, legend1, empty_plot2, legend2, empty_plot3, 
                     ncol = 1, 
                     align = 'hv',
                     rel_heights = c(0.2, 1, 0.3, 1, 0.6))  # Adjust these values to fine-tune the positions

# Arrange plots and combined legend
class_ps_class_plots <- plot_grid(plt_dendr_class_ps_class, plt_rel_class_ps_class, legends, 
                                  align = 'hv', 
                                  rel_widths = c(0.5, 1.5, 0.5),
                                  ncol = 3)
class_ps_class_plots


# This requires a directory called figures to be in your working directory
# All 3 lines have to be run for the file to be created.
png("figures/rarefied_Fig_Shotgun_vs_TE_relative_abundance_bySampleGroup.png", width = 2400, height = 1400, res = 100)
class_ps_class_plots
dev.off()
