# Canine RA plot

# =============================================================================
# CANINE FECAL ONLY - Relative abundance plot with dendrogram
# =============================================================================

# Subset to only Canine fecal samples
canine_group.ps <- subset_samples(
  all_AMR_group.ps,
  SampleGroup == "Canine fecal")

canine_group.ps <- prune_taxa(taxa_sums(canine_group.ps) > 0, canine_group.ps)

# Update InputDNA palette for canine (no Shotgun)
dna_pal_canine <- c(`50ng` = "#ffd700", `200ng` = "#e26700")  # adjust levels as needed

#############################################################################################
####################################   CSS TRANSFORM    ###################################

canine_group.ps.css <- phyloseq_transform_css(canine_group.ps, log = F)
canine_group.ps.css.df <- as(sample_data(canine_group.ps.css), "data.frame")

#############################################################################################
##############################         RELA ABUNDANCE         ###############################
#############################################################################################

rel_abund_canine <- transform_sample_counts(canine_group.ps.css, function(x) {x/sum(x)}*100)

# agglomerate
ra_class_canine <- tax_glom(rel_abund_canine, taxrank = "class")
ra_class_canine

ra_class_melt_temp_canine <- psmelt(ra_class_canine)

######## Re-label low abundance taxa into single category

# Convert 'class' to a character vector
ra_class_melt_temp_canine$class <- as.character(ra_class_melt_temp_canine$class)

# Step 1: Calculate median abundance for each class across all samples
all_sample_factor_by_abund_canine <- ra_class_melt_temp_canine %>%
  dplyr::group_by(class) %>%
  dplyr::summarize(median_class = median(Abundance)) %>%
  arrange(-median_class)

# Step 2: Identify classes where the median relative abundance is <= 1%
remainder_canine <- all_sample_factor_by_abund_canine %>%
  dplyr::filter(median_class <= 0.5) %>%
  dplyr::pull(class)

# Step 3: Relabel low-abundance classes as 'Low abundance classes (<1%)'
ra_class_melt_temp_canine$class <- ifelse(ra_class_melt_temp_canine$class %in% remainder_canine, 
                                          'Low abundance classes (<0.5%)', 
                                          ra_class_melt_temp_canine$class)

# Step 4: Group by Sample and relabeled 'class', sum Abundance
ra_class_melt_temp_canine <- ra_class_melt_temp_canine %>%
  dplyr::group_by(Sample, class) %>%
  dplyr::summarize(Abundance = sum(Abundance), .groups = "drop")

length(unique(ra_class_melt_temp_canine$class))

##
### Clustering at ASV level ####
##

# Cluster using hclust()
ps_AMR_canine.dist <- vegdist(t(otu_table(canine_group.ps.css)), method = "bray")
ps_AMR_canine.hclust <- hclust(ps_AMR_canine.dist)
plot(ps_AMR_canine.hclust) # example plot

# Extract data as dendrogram
ps_AMR_canine.dendro <- as.dendrogram(ps_AMR_canine.hclust)
ps_AMR_canine.dendro.data <- dendro_data(ps_AMR_canine.dendro, type = "rectangle")

# Sample names in order based on clustering
sample_names_ps_AMR_canine <- ps_AMR_canine.dendro.data$labels$label

# Add metadata
ps_AMR_canine.dendro.data$labels <- left_join(ps_AMR_canine.dendro.data$labels, 
                                              canine_group.ps.css.df, 
                                              by = c("label" = "Sample"))

# Setup the data, so that the layout is inverted
segment_data_ps_AMR_canine <- with(
  segment(ps_AMR_canine.dendro.data),
  data.frame(x = y, y = x, xend = yend, yend = xend))

# Use the dendrogram label data to position the gene labels
gene_pos_table_ps_AMR_canine <- with(
  ps_AMR_canine.dendro.data$labels, 
  data.frame(y_center = x, gene = as.character(label), x = y,
             InputDNA = as.character(InputDNA),
             height = 1))

gene_pos_table_ps_AMR_canine$InputDNA <- factor(gene_pos_table_ps_AMR_canine$InputDNA, 
                                                levels = c("50ng", "200ng"))  # adjust levels as needed

# Table to position the samples
sample_pos_table_ps_AMR_canine <- data.frame(sample = sample_names_ps_AMR_canine) %>%
  dplyr::mutate(x_center = (1:dplyr::n()), width = 1)

##
######## Relative abundance bar plot #########
##

ra_class_melt_canine <- ra_class_melt_temp_canine

# Use class melted data and add gene locations
joined_ra_class_ps_class_melt_canine <- ra_class_melt_canine %>%
  left_join(gene_pos_table_ps_AMR_canine, by = c("Sample" = "gene")) %>%
  left_join(sample_pos_table_ps_AMR_canine, by = c("Sample" = "sample"))

# Limits for the vertical axes
gene_axis_limits_ps_class_canine <- with(
  gene_pos_table_ps_AMR_canine, 
  c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))) + 0.1 * c(-1, 1)

# =============================================================================
# USE THE SAME COLOR PALETTE FROM THE MAIN FIGURE
# =============================================================================

# Get unique classes in canine data
canine_class_names <- unique(ra_class_melt_canine$class)

# Use the same named_palette from the main figure
# Make sure any classes not in the original palette get assigned colors
# First, check which classes are in canine but not in original palette
missing_classes <- setdiff(canine_class_names, names(named_palette))

if(length(missing_classes) > 0) {
  # Assign colors to missing classes from the remaining palette
  n_existing <- length(named_palette)
  additional_colors <- setNames(my_palette[(n_existing + 1):(n_existing + length(missing_classes))], 
                                missing_classes)
  named_palette_canine <- c(named_palette, additional_colors)
} else {
  named_palette_canine <- named_palette
}

# Make sure "Low abundance classes (<1%)" is gray
named_palette_canine["Low abundance classes (<0.5%)"] <- "gray50"

# Reorder class factors to match the main figure order where possible
# Get the order from the main figure
main_class_order <- names(named_palette)

# Classes present in canine data, in the same order as main figure
canine_classes_ordered <- main_class_order[main_class_order %in% canine_class_names]

# Add any classes unique to canine (shouldn't be many)
canine_only_classes <- setdiff(canine_class_names, main_class_order)
canine_classes_ordered <- c(canine_classes_ordered, canine_only_classes)

# Make sure "Low abundance classes (<1%)" is last
canine_classes_ordered <- c(setdiff(canine_classes_ordered, "Low abundance classes (<0.5%)"), 
                            "Low abundance classes (<0.5%)")

joined_ra_class_ps_class_melt_canine$class <- factor(joined_ra_class_ps_class_melt_canine$class,
                                                     levels = canine_classes_ordered,
                                                     ordered = TRUE)

#
## Relative abundance plot - Canine
#
plt_rel_class_ps_class_canine <- ggplot(joined_ra_class_ps_class_melt_canine, 
                                        aes(x = x_center, y = Abundance, fill = class, 
                                            height = height, width = width)) + 
  coord_flip() +
  geom_bar(stat = "identity", colour = "black", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = named_palette_canine, name = "AMR class", drop = FALSE) +
  scale_x_continuous(breaks = sample_pos_table_ps_AMR_canine$x_center, 
                     labels = sample_pos_table_ps_AMR_canine$sample, 
                     expand = c(0, 0)) + 
  labs(x = "", y = "Relative abundance") +
  theme_bw() +
  theme(legend.position = "right",
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.75),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.75, lineend = "square", colour = "black"),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, colour = "black"),
        plot.margin = unit(c(1, 0.2, 0.2, -0.7), "cm"), 
        panel.grid.minor = element_blank())

plt_rel_class_ps_class_canine

## Dendrogram plot - Canine (only InputDNA annotation) ####

bar_w <- 0.04
gene_pos_table_ps_AMR_canine <- gene_pos_table_ps_AMR_canine %>%
  dplyr::mutate(x_dna = x)

plt_dendr_class_ps_class_canine <- ggplot(segment_data_ps_AMR_canine) +
  ## dendrogram branches
  geom_segment(aes(x, y, xend = xend, yend = yend),
               linewidth = .75, lineend = "round", linejoin = "round") +
  ## Only InputDNA tiles
  geom_tile(data = gene_pos_table_ps_AMR_canine,
            aes(x = x_dna, y = y_center, fill = InputDNA),
            width = bar_w, height = 1, colour = NA) +
  scale_fill_manual("Input DNA (ng)", values = dna_pal_canine, drop = FALSE) +
  scale_y_continuous(breaks = gene_pos_table_ps_AMR_canine$y_center,
                     labels = gene_pos_table_ps_AMR_canine$gene,
                     limits = gene_axis_limits_ps_class_canine,
                     expand = c(0, 0)) +
  scale_x_reverse() +
  labs(x = "Bray-Curtis distance", y = NULL) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = "black", linewidth = 0.75),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        legend.position = "right")

print(plt_dendr_class_ps_class_canine)

## Combine relative abundance plot and dendrogram - Canine (FIXED ALIGNMENT)

# Extract legends
legend1_canine <- cowplot::get_legend(plt_dendr_class_ps_class_canine)
legend2_canine <- cowplot::get_legend(plt_rel_class_ps_class_canine)

# Remove legends from the original plots AND fix margins
plt_dendr_class_ps_class_canine_fixed <- plt_dendr_class_ps_class_canine + 
  theme(legend.position = "none",
        plot.margin = unit(c(1, 0, 0.2, 0.2), "cm"))

plt_rel_class_ps_class_canine_fixed <- plt_rel_class_ps_class_canine + 
  theme(legend.position = "none",
        plot.margin = unit(c(1, 0.2, 0.2, 0), "cm"))

# Make sure both plots use the same y-axis limits
plt_dendr_class_ps_class_canine_fixed <- plt_dendr_class_ps_class_canine_fixed +
  scale_y_continuous(breaks = gene_pos_table_ps_AMR_canine$y_center,
                     labels = gene_pos_table_ps_AMR_canine$gene,
                     limits = gene_axis_limits_ps_class_canine,
                     expand = c(0, 0))

plt_rel_class_ps_class_canine_fixed <- plt_rel_class_ps_class_canine_fixed +
  scale_x_continuous(breaks = sample_pos_table_ps_AMR_canine$x_center,
                     labels = sample_pos_table_ps_AMR_canine$sample,
                     limits = gene_axis_limits_ps_class_canine,  # Use same limits!
                     expand = c(0, 0))

## make empty plots to help with spacing legends
empty_plot1 <- ggplot() + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
empty_plot2 <- ggplot() + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
empty_plot3 <- ggplot() + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))

# Arrange legends
legends_canine <- plot_grid(empty_plot1, legend1_canine, empty_plot2, legend2_canine, empty_plot3, 
                            ncol = 1, 
                            rel_heights = c(0.2, 0.5, 0.3, 1, 0.6))

# Arrange plots with proper alignment
class_ps_class_plots_canine <- plot_grid(
  plt_dendr_class_ps_class_canine_fixed, 
  plt_rel_class_ps_class_canine_fixed, 
  legends_canine, 
  align = 'h',
  axis = 'tb',
  rel_widths = c(0.3, 1, 0.5),
  ncol = 3
)

class_ps_class_plots_canine

# Save figure
png("figures/Fig_Canine_relative_abundance_byInputDNA.png", width = 1700, height = 800, res = 100)
class_ps_class_plots_canine
dev.off()
