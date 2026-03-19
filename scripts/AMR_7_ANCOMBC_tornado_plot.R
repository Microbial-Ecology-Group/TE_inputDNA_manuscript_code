library(ANCOMBC)
library(dplyr)
library(DT)
library(RColorBrewer)

DatasetA.group_AMR.ps <- subset_samples(group.ps, SampleGroup != "Canine fecal")

colnames(phyloseq::tax_table(DatasetA.group_AMR.ps))[4] <- "Species"

output_Species_datasetA = ancombc2(data = DatasetA.group_AMR.ps, tax_level = "class",
                             fix_formula ="Kit + InputDNA + OriginalSample", rand_formula = NULL,
                             p_adj_method = "holm", 
                             prv_cut = 0.25, lib_cut = 1000, s0_perc = 0.05,
                             group = "Kit", struc_zero = TRUE, neg_lb = TRUE,
                             alpha = 0.05, n_cl = 2, verbose = TRUE,
                             global = TRUE, pairwise = TRUE, dunnet = FALSE, trend = FALSE)



# Extract the pairwise test results
pairwise_ancom_output_DatasetA <- output_Species_datasetA$res

# Check the structure of the pairwise results
str(pairwise_ancom_output_DatasetA)

# Filter the results to focus on comparisons between Kits
# For example, if you want to look at a specific species or taxon, you can filter accordingly
df_pairwise_ancom_output_DatasetA = pairwise_ancom_output_DatasetA %>%
  dplyr::select(taxon, contains("Kit")) 

df_fig_Treatment_DatasetA <- df_pairwise_ancom_output_DatasetA %>%
  # Step 1: Apply the condition for the lfc_Kit columns
  dplyr::mutate(across(starts_with("lfc_Kit"), ~ifelse(get(gsub("lfc_", "diff_", cur_column())) == 1, .x,.x)),
                # Step 2: Round the lfc_Kit columns and create _rounded columns
                across(starts_with("lfc_Kit"), ~round(.x, 3), .names = "{.col}_rounded"),
                # Step 3: Assign colors based on diff_Kit and passed_ss columns
                across(starts_with("diff_Kit"), 
                       ~ifelse(.x == 1 & get(gsub("diff_", "passed_ss_", cur_column())) == 1, "darkred", "darkgrey"), 
                       .names = "{.col}_color")) %>%
  # Step 4: Pivot the _rounded columns into long format
  pivot_longer(cols = contains("_rounded"), names_to = "group", values_to = "value", names_prefix = "lfc_") %>%
  # Step 5: Pivot the _color columns into long format
  pivot_longer(cols = contains("_color"), names_to = "color_group", values_to = "color", names_prefix = "diff_") %>%
  # Step 6: Filter to ensure the group and color_group match
  filter(gsub("_rounded", "", group) == gsub("_color", "", color_group)) %>%
  # Step 7: Select the relevant columns and arrange by taxon
  select(-color_group) %>%
  arrange(taxon)


df_fig_Treatment_DatasetA_filtered <- df_fig_Treatment_DatasetA %>%
  dplyr::filter(str_ends(group, "_rounded")) %>%
  dplyr::mutate(group = dplyr::case_when(
    group == "KitXT HS_rounded" ~ "XT HS vs XT",
    TRUE ~ group)) %>%
  dplyr::filter(group %in% c("XT HS vs XT")) %>%
  dplyr::mutate(group = factor(group, levels = c("XT HS vs XT"))) %>%
  droplevels()




levels(as.factor(df_fig_Treatment_DatasetA_filtered$group))

#### DatasetA Dotplot - but with bias corrected abundances #######

# Assuming 'ancom_output_DatasetA$bias_correct_log_table' is your bias-corrected log-transformed data
log_table_DatasetA <- output_Species_datasetA$bias_correct_log_table

# Replace NA values in the log-transformed data with 0
log_table_DatasetA[is.na(log_table_DatasetA)] <- 0

# Step 1: Exponentiate the log-transformed data to get values on the original scale
pseudo_counts_DatasetA <- exp(log_table_DatasetA)

# Step 2: Scale the pseudo-counts (e.g., normalize each sample to a total count)
# For example, you can scale the pseudo-counts to a total library size or another baseline.

# Calculate the total sum of counts for each sample before scaling
total_counts_DatasetA <- apply(pseudo_counts_DatasetA, 2, sum)

# Define a library size (or use the sum of pseudo_counts as the new total count)
desired_total_count_DatasetA <- 100000  # For example, set each sample to have a total of 10,000 counts

# Scale the pseudo-counts to the desired total count per sample
scaled_counts_DatasetA <- t(t(pseudo_counts_DatasetA) * (desired_total_count_DatasetA / total_counts_DatasetA))

# Convert scaled_counts_DatasetA to a data frame
scaled_counts_DatasetA.df <- as.data.frame(scaled_counts_DatasetA)

# Add a column for taxa names
scaled_counts_DatasetA.df$taxon <- rownames(scaled_counts_DatasetA.df)

# Convert from wide to long format using pivot_longer
scaled_counts_DatasetA.long <- scaled_counts_DatasetA.df %>%
  pivot_longer(cols = -taxon, names_to = "Sample", values_to = "BiasAdj_counts")

# Make melted data
DatasetA.AMR_data_group.melt <- psmelt(DatasetA.group_AMR.ps)

# Merge data
# Merge the melted data with the scaled counts on 'taxon' and 'Sample'
BiasAdj_DatasetA.AMR_data_group.melt <- DatasetA.AMR_data_group.melt %>%
  left_join(scaled_counts_DatasetA.long, by = c("class" = "taxon", "Sample" = "Sample"))

Species_avg_counts_DatasetA_by_Kit <- BiasAdj_DatasetA.AMR_data_group.melt %>%
  dplyr::group_by(class, Kit) %>%
  dplyr::summarise(Avg_BiasAdj_Count = mean(BiasAdj_counts, na.rm = FALSE)) %>%
  ungroup() %>%
  dplyr::mutate(Kit_label = dplyr::case_when(
    Kit == 'XT HS' ~ "XT HS vs XT",
    TRUE ~ as.character(Kit)  # Keep other Kit values unchanged
  ))


df_fig_Treatment_DatasetA_filtered_merged_data <- df_fig_Treatment_DatasetA_filtered %>%
  left_join(Species_avg_counts_DatasetA_by_Kit, by = c("group" = "Kit_label", "taxon" = "class"))

# Extract taxonomic data from the phyloseq object
tax_data_DatasetA <- as.data.frame(tax_table(DatasetA.group_AMR.ps))

# Join the taxonomic data with your data frame
df_with_taxonomy_DatasetA <- df_fig_Treatment_DatasetA_filtered_merged_data %>%
  left_join(tax_data_DatasetA, by = c("taxon" = "class"),relationship = "many-to-many")  # Assuming 'Species' is your key column to join

# # Reorder the taxon based on class, mechanism, and species
# df_with_taxonomy_DatasetA <- df_with_taxonomy_DatasetA %>%
#   arrange(class, mechanism, taxon) %>%  # Sort the data by class, mechanism, and Species
#   dplyr::mutate(taxon = factor(taxon, levels = rev(unique(taxon))))  # Reorder the 'taxon' factor based on the new order

df_with_taxonomy_DatasetA <- df_with_taxonomy_DatasetA %>%
  # Arrange the data by value (for ordering taxon based on the value column)
  arrange(value) %>%
  # Reorder the taxon factor based on the values in the value column
  dplyr::mutate(taxon = factor(taxon, levels = unique(taxon[order(-value)])))

#df_with_taxonomy_DatasetA$group = factor(df_with_taxonomy_DatasetA$group, levels = c("From Days 0&1 to Days 5&6", "From Days 5&6 to Days 14&15","From Days 14&15 to Days 20&21"))


## Trying to add taxonomy
main_species_plot_DatasetA <- ggplot(df_with_taxonomy_DatasetA, aes(x = value, y = taxon, color = color)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "solid") +  # Add vertical dashed line at 0
  geom_point(aes(size = Avg_BiasAdj_Count)) +  # Set the size based on Avg_BiasAdj_Count
  scale_color_identity() +  # Use color directly from the dataset
  scale_size_continuous(range = c(1, 6)) +  # Increase the size range to create more contrast
  labs(x = "Log-Fold Change", y = "Taxon") +
  theme_minimal() +  # Set the theme to minimal
  theme(
    axis.text.y =element_text(size = 10, colour = "black"),  # Hide y-axis species labels (to show in taxonomy plot)
    axis.text.x = element_text(size = 14, angle = 45, vjust = 0.5),
    axis.title.y = element_blank(),  # Remove y-axis title
    axis.title.x = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),  # Remove strip background
    strip.text = element_text(size = 12, colour = "black"), 
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Centered title
    legend.position = "none",  # Hide the legend as the color is used directly
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),  # Make panel background transparent
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0),  # Add border around facets
    panel.spacing = unit(1, "lines")  # Add space between facets
  ) +
  facet_wrap(~group, ncol = 3, scales = "fixed")


main_species_plot_DatasetA


# Save the final plot
ggsave("../Writing/Final_figures/Figure_7_DA_DatasetA_by_Kit.jpg", plot = main_species_plot_DatasetA, width = 10, height = 11, dpi = 300)








# Create the taxonomy plot data
taxonomy_plot_data_DatasetA <- df_with_taxonomy_DatasetA %>%
  distinct(class, mechanism, taxon)

# Modify the data to create new columns with the "label_" prefix
taxonomy_plot_data_DatasetA <- taxonomy_plot_data_DatasetA %>%
  dplyr::group_by(class) %>%
  dplyr::mutate(label_class = ifelse(row_number() == 1, class, "")) %>%  # Create 'label_class' with only the first occurrence of each class
  ungroup() %>%
  dplyr::group_by(class, mechanism) %>%
  dplyr::mutate(label_mechanism = ifelse(row_number() == 1, mechanism, ""))  # Create 'label_mechanism' with only the first occurrence of each mechanism

# Create the updated taxonomy plot
taxonomy_plot_DatasetA <- ggplot(taxonomy_plot_data_DatasetA) +
  geom_text(aes(x = 1, y = taxon, label = label_class), hjust = 1, size = 3.5) +  # Use 'label_class' for class
  geom_text(aes(x = 2, y = taxon, label = label_mechanism), hjust = 1, size = 3.5) +  # Use 'label_mechanism' for mechanism
  scale_x_continuous(limits = c(0.5, 2.5), breaks = c(1, 2), labels = c("Class", "Mechanism")) +
  scale_y_discrete(limits = levels(taxonomy_plot_data_DatasetA$taxon)) +  # Ensure y-axis matches the taxon order
  theme_void() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")  # Reduce margin on the right of taxonomy_plot
  )


combined_plot_DatasetA  <- plot_grid(
  taxonomy_plot_DatasetA  + theme(plot.margin = unit(c(1, 0, 0, 2), "cm")),  # Remove margins from the taxonomy plot
  main_species_plot_DatasetA  + theme(plot.margin = unit(c(1, 2, 0, 0), "cm")),  # Remove margins from the main plot
  ncol = 2, 
  rel_widths = c(1.9,1.8),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb"
)


# Save the final plot
ggsave("figures/DA_DatasetA_by_Kit.jpg", plot = combined_plot_DatasetA, width = 20, height = 11, dpi = 300)




