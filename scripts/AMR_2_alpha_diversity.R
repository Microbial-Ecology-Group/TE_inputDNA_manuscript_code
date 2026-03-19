#############################################################################################
##############################         ALPHA DIVERSITY         ##############################
#############################################################################################
#############################################################################################


alpha_div <- estimate_richness(group.ps, measures = c("Observed","Shannon","Simpson","InvSimpson"))
alpha_div
alpha_div.df <- as(sample_data(group.ps), "data.frame")
alpha_div_meta <- cbind(alpha_div, alpha_div.df)

write.csv(alpha_div_meta, "agilent_alpha_div_measures.csv")

## Statistical tests ####
# Kit InputDNA      SampleGroup Dilution

library(purrr)



# Example usage
# Assuming your data frame is 'alpha_div_meta' and you want to run the test on 'Observed' for each factor in 'SampleGroup'

# Run the nested pairwise Wilcoxon test function
nested_Observed_Kit_results <- nested_pairwise_wilcox(
  data = alpha_div_meta,
  measure_col = "Observed",    # The column with the values you're testing
  group_col = "Kit",   # The column defining the groups for comparison
  nesting_factor = "SampleGroup" # The column by which to nest the tests
)

# Inspect the Kit_results (each element in the list corresponds to a test result for a SampleGroup factor)
nested_Observed_Kit_results


# Run the nested pairwise Wilcoxon test function
nested_Shannon_Kit_results <- nested_pairwise_wilcox(
  data = alpha_div_meta,
  measure_col = "Shannon",    # The column with the values you're testing
  group_col = "Kit",   # The column defining the groups for comparison
  nesting_factor = "SampleGroup" # The column by which to nest the tests
)

# Inspect the Kit_results (each element in the list corresponds to a test result for a SampleGroup factor)
nested_Shannon_Kit_results

# Run the nested pairwise Wilcoxon test function
nested_Observed_InputDNA_results <- nested_pairwise_wilcox(
  data = alpha_div_meta,
  measure_col = "Observed",    # The column with the values you're testing
  group_col = "InputDNA",   # The column defining the groups for comparison
  nesting_factor = "SampleGroup" # The column by which to nest the tests
)

# Inspect the InputDNA_results (each element in the list corresponds to a test result for a SampleGroup factor)
nested_Observed_InputDNA_results


# Run the nested pairwise Wilcoxon test function
nested_Shannon_InputDNA_results <- nested_pairwise_wilcox(
  data = alpha_div_meta,
  measure_col = "Shannon",    # The column with the values you're testing
  group_col = "InputDNA",   # The column defining the groups for comparison
  nesting_factor = "SampleGroup" # The column by which to nest the tests
)

# Inspect the InputDNA_results (each element in the list corresponds to a test result for a SampleGroup factor)
nested_Shannon_InputDNA_results


## Dataset A ####
DatasetA_metadata <- subset(alpha_div_meta, SampleGroup != "Canine fecal")

## Observed summary
mean(DatasetA_metadata$Observed)
min(DatasetA_metadata$Observed)
max(DatasetA_metadata$Observed)

## Shannon summary
mean(DatasetA_metadata$Shannon)
min(DatasetA_metadata$Shannon)
max(DatasetA_metadata$Shannon)


### SampleGroup
SampleGroup.richness.pw <- pairwise.wilcox.test(DatasetA_metadata$Observed, DatasetA_metadata$SampleGroup, p.adjust.method = "BH")
SampleGroup.richness.pw # p-values in the matrix

SampleGroup.shannon.pw <- pairwise.wilcox.test(DatasetA_metadata$Shannon, DatasetA_metadata$SampleGroup, p.adjust.method = "BH")
SampleGroup.shannon.pw # p-values in the matrix

# Library type
Kit.richness.pw <- pairwise.wilcox.test(DatasetA_metadata$Observed, DatasetA_metadata$Kit, p.adjust.method = "BH")
Kit.richness.pw # p-values in the matrix

Kit.shannon.pw <- pairwise.wilcox.test(DatasetA_metadata$Shannon, DatasetA_metadata$Kit, p.adjust.method = "BH")
Kit.shannon.pw # p-values in the matrix


### Kit
Kit.richness.pw <- pairwise.wilcox.test(DatasetA_metadata$Observed, DatasetA_metadata$Kit, p.adjust.method = "BH")
Kit.richness.pw # p-values in the matrix

Kit.shannon.pw <- pairwise.wilcox.test(DatasetA_metadata$Shannon, DatasetA_metadata$Kit, p.adjust.method = "BH")
Kit.shannon.pw # p-values in the matrix


### By kit #####

XT_metadata <- subset(alpha_div_meta, Kit == "XT")

InputDNA.richness.pw <- pairwise.wilcox.test(XT_metadata$Observed, XT_metadata$InputDNA, p.adjust.method = "BH")
InputDNA.richness.pw # p-values in the matrix

InputDNA.shannon.pw <- pairwise.wilcox.test(XT_metadata$Shannon, XT_metadata$InputDNA, p.adjust.method = "BH")
InputDNA.shannon.pw # p-values in the matrix

XT_HS_metadata <- subset(alpha_div_meta, Kit == "XT_HS")

InputDNA.richness.pw <- pairwise.wilcox.test(XT_HS_metadata$Observed, XT_HS_metadata$InputDNA, p.adjust.method = "BH")
InputDNA.richness.pw # p-values in the matrix

InputDNA.shannon.pw <- pairwise.wilcox.test(XT_HS_metadata$Shannon, XT_HS_metadata$InputDNA, p.adjust.method = "BH")
InputDNA.shannon.pw # p-values in the matrix



### Dataset B #####

DatasetB_metadata <- subset(alpha_div_meta, SampleGroup == "Canine fecal")

# Input DNA

InputDNA.richness.pw <- pairwise.wilcox.test(DatasetB_metadata$Observed, DatasetB_metadata$InputDNA, p.adjust.method = "BH")
InputDNA.richness.pw # p-values in the matrix

InputDNA.shannon.pw <- pairwise.wilcox.test(DatasetB_metadata$Shannon, DatasetB_metadata$InputDNA, p.adjust.method = "BH")
InputDNA.shannon.pw # p-values in the matrix


## Dilution Figures ####

ggplot(alpha_div_meta, aes(x= Dilution, y= Observed, fill = Dilution, colour = Dilution)) +
  theme_bw() + 
  labs(y= "Observed ARGs") +
  geom_boxplot(alpha=0.4) +
  geom_point() +
  theme(legend.position = "right",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size =24, colour = "white"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 28),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 28),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())


dilution.richness.pw <- pairwise.wilcox.test(alpha_div_meta$Observed, alpha_div_meta$Dilution, p.adjust.method = "BH")
dilution.richness.pw # p-values in the matrix



ggplot(alpha_div_meta, aes(x= Dilution, y= Shannon, fill = Dilution, colour = Dilution)) +
  theme_bw() + 
  labs(y= "Shannon's diversity") +
  geom_boxplot(alpha=0.4) +
  geom_point() +
  theme(legend.position = "right",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size =24, colour = "white"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 28),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 28),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

dilution.shannon.pw <- pairwise.wilcox.test(alpha_div_meta$Shannon, alpha_div_meta$Dilution, p.adjust.method = "BH")
dilution.shannon.pw # p-values in the matrix



### Figure 3 SampleGroup - Observed ####

plot_div_obs <- ggplot(alpha_div_meta[alpha_div_meta$SampleGroup != "Canine fecal",], aes(x= SampleGroup, y= Observed, fill = Kit, colour = Kit)) +
  theme_bw() + 
  labs(y= "Observed AMR groups") +
  geom_boxplot(alpha=0.4) +
  geom_point(position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = c( "darkgreen", "maroon4")) +
  scale_color_manual(values = c( "darkgreen", "maroon4")) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size =24, colour = "white"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 28),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 28),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())


SampleGroup.richness.pw <- pairwise.wilcox.test(alpha_div_meta$Observed, alpha_div_meta$SampleGroup, p.adjust.method = "BH")
SampleGroup.richness.pw # p-values in the matrix



plot_div_shan <- ggplot(alpha_div_meta[alpha_div_meta$SampleGroup != "Canine fecal",], aes(x= SampleGroup, y= Shannon, fill = Kit, colour = Kit)) +
  theme_bw() + 
  labs(y= "Shannon's diversity") +
  geom_boxplot(alpha=0.4) +
  geom_point(position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = c( "darkgreen", "maroon4")) +
  scale_color_manual(values = c("darkgreen", "maroon4")) +
  theme(legend.position = "right",
        legend.text = element_text(size = 18),    # Increase legend text size
        legend.title = element_text(size = 20),   # Increase legend title size
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size = 24, colour = "white"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 28),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 28),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())


SampleGroup.shannon.pw <- pairwise.wilcox.test(alpha_div_meta$Shannon, alpha_div_meta$SampleGroup, p.adjust.method = "BH")
SampleGroup.shannon.pw # p-values in the matrix


diversity_grid <- plot_grid(
  plot_div_obs,
  plot_div_shan,
  ncol = 2,
  labels = c("A", "B")
)


# Save the final plot
ggsave("../Writing/Final_figures/Fig3_diversity_boxplot.jpg", plot = diversity_grid, width = 20, height = 11, dpi = 300)


# Kit
pairwise.wilcox.test(alpha_div_meta$Observed, alpha_div_meta$Kit, p.adjust.method = "BH")
pairwise.wilcox.test(alpha_div_meta$Shannon, alpha_div_meta$Kit, p.adjust.method = "BH")

# Dilution
pairwise.wilcox.test(alpha_div_meta$Observed, alpha_div_meta$Dilution, p.adjust.method = "BH")
pairwise.wilcox.test(alpha_div_meta$Shannon, alpha_div_meta$Dilution, p.adjust.method = "BH")

# InputDNA
pairwise.wilcox.test(alpha_div_meta$Observed, alpha_div_meta$InputDNA, p.adjust.method = "BH")
pairwise.wilcox.test(alpha_div_meta$Shannon, alpha_div_meta$InputDNA, p.adjust.method = "BH")

mean_values <- alpha_div_meta %>%
  group_by(InputDNA) %>%
  summarise(Mean_Shannon = mean(Shannon, na.rm = TRUE))

mean(alpha_div_meta$Shannon)
min(alpha_div_meta$Shannon)
max(alpha_div_meta$Shannon)


boxplot(alpha_div_meta$Shannon ~alpha_div_meta$InputDNA)


# SampleGroup
pairwise.wilcox.test(alpha_div_meta$Observed, alpha_div_meta$SampleGroup, p.adjust.method = "BH")
pairwise.wilcox.test(alpha_div_meta$Shannon, alpha_div_meta$SampleGroup, p.adjust.method = "BH")

boxplot(alpha_div_meta$Shannon ~alpha_div_meta$SampleGroup)

## Figure - by InputDNA shape by Kit ####
