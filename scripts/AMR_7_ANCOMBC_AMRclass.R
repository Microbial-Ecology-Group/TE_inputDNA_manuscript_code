# ANCOMBC - class level
library(ANCOMBC)
library(dplyr)
library(DT)
library(RColorBrewer)

DatasetA.data_AMR.ps <- subset_samples(data_AMR.ps, SampleGroup != "Canine fecal")
DatasetA.data_AMR.class.ps <- tax_glom(DatasetA.data_AMR.ps, taxrank = "class")

colnames(phyloseq::tax_table(DatasetA.data_AMR.class.ps))[2] <- "Kingdom"

colnames(phyloseq::tax_table(DatasetA.data_AMR.ps))[2] <- "Kingdom"

DatasetA.data_AMR.ps
output_group_datasetA = ancombc2(data = DatasetA.data_AMR.ps, tax_level = "Kingdom",
                                 fix_formula ="LIbraryType + InputDNA  + OriginalSample", rand_formula = NULL,
                                 p_adj_method = "holm", 
                                 prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                 group = "LIbraryType", struc_zero = TRUE, neg_lb = TRUE,
                                 alpha = 0.05, n_cl = 2, verbose = TRUE,
                                 global = TRUE, pairwise = TRUE, dunnet = FALSE, trend = FALSE)


## Structural zeros ####
tab_zero = output_group_datasetA$zero_ind
tab_zero %>%
  datatable(caption = "The detection of structural zeros")

## Sensitivity scores ####
tab_sens = output_group_datasetA$pseudo_sens_tab
tab_sens %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(tab_sens)[-1], digits = 2)



### ANCOM-BC2 primary analysis ####
#Result from the ANCOM-BC2 methodology to determine taxa that are differentially abundant
#according to the covariate of interest. Results contain: 1) log fold changes, 2) standard errors, 
#3) test statistics, 4) p-values, 5) adjusted p-values, 6) indicators of whether the taxon is
#differentially abundant (TRUE) or not (FALSE).

res_prim = output_group_datasetA$res

## Results for LIbraryType order  #####
# Original results
df_LIbraryType = res_prim %>%
  dplyr::select(taxon, contains("LIbrary")) 


df_fig_LIbraryType = df_LIbraryType %>%
  filter(diff_LIbraryTypeXT_HS == 1) %>%
  mutate(lfc_LIbraryType = ifelse(diff_LIbraryTypeXT_HS == 1, 
                                  lfc_LIbraryTypeXT_HS, 0)) %>%
  transmute(taxon, `XT_HS vs XT` = round(lfc_LIbraryType, 2)) %>%
  pivot_longer(cols = `XT_HS vs XT`, names_to = "group", values_to = "value") %>%
  arrange(-value) %>%
  mutate(qvalue = df_LIbraryType$q_LIbraryTypeXT_HS[match(taxon, df_LIbraryType$taxon)])


# Only one direction of values
# fig_LIbraryType = df_fig_LIbraryType %>%
#   ggplot(aes(x = group, y = reorder(taxon, value), fill = value)) + 
#   geom_tile(color = "black") +
#   scale_fill_gradientn(colors =  rev(brewer.pal(9, "Blues")[c(1:8)]),
#                        values = rescale(-c(up, lo), to = c(0, 1)),
#                        na.value = "white", 
#                        name = NULL) +
#   geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
#   labs(x = NULL, y = NULL, title = "LogFC in dogs compared by LIbraryType type") +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5))
# fig_LIbraryType

fig_ANCOMBC_LIbraryType = df_fig_LIbraryType %>%
  ggplot(aes(x = group, y = reorder(taxon, value), fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")),
                       limits = range(df_fig_LIbraryType$value),
                       na.value = "white", 
                       name = "logFC") +
  geom_text(aes(group, taxon, label = value), color = "white", size = 4) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "left",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14) )

fig_ANCOMBC_LIbraryType

### Taxa plots ####

# Extract the results data frame from the ANCOMBC2 output
results_df = output_group_datasetA$res

# Filter the results to only include significant features
significant_features = df_fig_LIbraryType$taxon

# Extract the names of the significant features
significant_feature_names = as.character(significant_features)

# Extract the bias_correct_log_table from the ancombc output
log_table = output_group_datasetA$bias_correct_log_table

# Subset the table to include only the significant features
log_table_significant = log_table[rownames(log_table) %in% significant_feature_names, ]

# Convert the row names into a column
log_table_significant$Feature = rownames(log_table_significant)

# Reshape the data into a format suitable for ggplot2
df = reshape2::melt(log_table_significant, id.vars = "Feature")

# Extract the sample data from the phyloseq object and convert it into a data frame
sample_data_df = data.frame(sample_data(DatasetA.data_AMR.class.ps))

# Add a column with the sample names
sample_data_df$Sample = rownames(sample_data_df)

# Merge the df data frame with the sample data
merged_df = merge(df, sample_data_df, by.x = "variable", by.y = "Sample")


# # Transform the phyloseq object into a data frame
# df = psmelt(DatasetA.data_AMR.class.ps)
# # Subset the data frame to only include the significant features
# df_significant = df[df$Kingdom %in% significant_feature_names, ]


merged_df$Class <- merged_df$Feature

merged_df$Class <- factor(merged_df$Class, levels = c("Rifampin","Elfamycins",
                                                                "Aminoglycosides", "Sodium resistance","Phenicol","Sulfonamides"))


# Create a boxplot of feature abundance by "Kingdom", compared by "LIbraryType"
# ggplot(df_significant, aes(x = Class, y = Abundance, fill = LIbraryType)) +
#   geom_boxplot(alpha = 0.6) +
#   geom_jitter(width = 0.35, size = .8) +
#   labs(x = "Class", y = "Abundance", fill = "LIbraryType") +
#   theme_minimal()

# Create the plot
fig_sigcounts_LIbraryType <- ggplot(merged_df, aes(x = Class, y = value, fill = LIbraryType)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.35, size = .8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_discrete(limits = levels(merged_df$Class)) +
  scale_fill_manual(values = c( "Gold","#008080")) +
  coord_flip() +
  labs(x = "Class", y = "Bias-corrected log abundance", fill = "LIbraryType") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_blank())

# Print the plot
print(fig_sigcounts_LIbraryType)


fig_ANCOMBC_grid <- plot_grid(fig_ANCOMBC_LIbraryType, fig_sigcounts_LIbraryType,labels = c("A", " B"), align = "h", ncol = 2,rel_widths = c(0.3, 0.7))
#title = "LogFC in samples processed with the XTHS vs the XT kit") 
##
### Structural zeros ####
##
tab_zero = output_group_datasetA$zero_ind
tab_zero %>%
  datatable(caption = "The detection of structural zeros")



# join the two dataframes by the taxon column
joined_df <- left_join(df_fig_LIbraryType, tab_zero, by = "taxon")

# check if any significantly different taxons have structural zeros
any(!is.na(joined_df$`structural_zero (LIbraryType = XT)`) | 
      !is.na(joined_df$`structural_zero (LIbraryType = XT_HS)`))



# Get taxons with structural zeros for XT LIbraryType type
structural_zero_XT <- tab_zero[tab_zero$`structural_zero (LIbraryType = XT)` == TRUE,]$taxon

# Get taxons with structural zeros for XT_HS LIbraryType type
structural_zero_XT_HS <- tab_zero[tab_zero$`structural_zero (LIbraryType = XT_HS)` == TRUE,]$taxon

# Check if significantly different taxons match the structural zeros
match_XT <- df_fig_LIbraryType$taxon %in% structural_zero_XT
match_XT_HS <- df_fig_LIbraryType$taxon %in% structural_zero_XT_HS

# Print taxons that have structural zeros in the ANCOMBC2 model
cat("Taxons with structural zeros for XT LIbraryType type:\n")
print(structural_zero_XT)
cat("\n")

cat("Taxons with structural zeros for XT_HS LIbraryType type:\n")
print(structural_zero_XT_HS)
cat("\n")

# Print whether the significantly different taxons match the structural zeros
cat("Do significantly different taxons match the structural zeros for XT LIbraryType type?\n")
print(all(match_XT == structural_zero_XT))

cat("Do significantly different taxons match the structural zeros for XT_HS LIbraryType type?\n")
print(all(match_XT_HS == structural_zero_XT_HS))





##
### Sensitivity scores ####
##

tab_sens = output_group_datasetA$ss_tab
tab_sens %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(tab_sens)[-1], digits = 2)


sens_type = tab_sens %>%
  transmute(taxon, sens_type = LIbraryTypeXT_HS) %>%
  left_join(df_LIbraryType %>%
              transmute(taxon, diff_type = diff_LIbraryTypeXT_HS), 
            by = "taxon") %>%
  mutate(group = "XT_HS vs. XT")

sens_type$diff_type = recode(sens_type$diff_type * 1, 
                             `1` = "Significant",
                             `0` = "Nonsignificant")

fig_sens_type = sens_type %>%
  ggplot(aes(x = taxon, y = sens_type, color = diff_type)) +
  geom_point() +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  facet_grid(rows = vars(group), scales = "free") +
  labs(x = NULL, y = "Sensitivity Score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5))
fig_sens_type



