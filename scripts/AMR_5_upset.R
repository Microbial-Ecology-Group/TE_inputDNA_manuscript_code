# Upset plots

#Taras Plots# Much simpler than Lee's but his still works

library(UpSetR)

library("MicrobiotaProcess")


group_AMR.ps <- tax_glom(data_AMR.ps, taxrank = "group")
class_AMR.ps <- tax_glom(data_AMR.ps, taxrank = "class")

## Apply prevalence filter ####
# Step 1: Melt the phyloseq object into a data frame
melted_df <- psmelt(group_AMR.ps)

# Step 2: Calculate prevalence for each SampleGroup
# Prevalence is the number of samples with a taxon > 0 divided by the total number of samples in that group
prevalence_df <- melted_df %>%
  dplyr::group_by(SampleGroup, OTU) %>%
  dplyr::summarise(prevalence = sum(Abundance > 0) / n()) %>%  # Calculate prevalence correctly
  ungroup()

# Step 3: Filter taxa that are present in at least 10% of samples in each SampleGroup
filtered_taxa <- prevalence_df %>%
  filter(prevalence >= 0.05) %>%
  pull(OTU)  # Get the names of the OTUs that meet the criteria

# Step 4: Filter the original phyloseq object based on the filtered taxa
filtered_phyloseq <- prune_taxa(filtered_taxa, group_AMR.ps)

# View the filtered phyloseq object
filtered_phyloseq


## Sample group ####

# samplegroup_upsetda_edit <- get_upset(filtered_phyloseq, factorNames="SampleGroup") ## ASV
# 
# 
# upset(samplegroup_upsetda_edit, sets=c("Avian fecal","Bovine fecal", "Human wastewater", "Porcine fecal", "Canine fecal"),
#       
#       sets.bar.color = c("#d6e2e9","#e9d556","#ac1d1c", "blue", "orange"),text.scale = 2,
#       
#       order.by = "freq", empty.intersections = "on")
# 

# Kit
libkit_unfiltered_upsetda_edit <- get_upset(group_AMR.ps, factorNames="Kit") ## ASV

group_level_upset_unfiltered <- upset(libkit_unfiltered_upsetda_edit, sets=c("XT","XT HS"),
                           
                           sets.bar.color = c("darkgreen", "maroon4"),text.scale = 2,
                           
                           order.by = "freq", empty.intersections = "on")

group_level_upset_unfiltered

# Kit
libkit_upsetda_edit <- get_upset(filtered_phyloseq, factorNames="Kit") ## ASV

group_level_upset <- upset(libkit_upsetda_edit, sets=c("XT","XT HS"),
      
      sets.bar.color = c("navyblue", "gold2"),text.scale = 2,
      
      order.by = "freq", empty.intersections = "on")

group_level_upset


# Then, subset the OTU table to include only rows corresponding to Pactamycin
Pactamycin_otus <- subset_taxa(class_AMR.ps, class == "Pactamycin")

# Now, convert the OTU table to a dataframe
Pactamycin_df <- as.data.frame(otu_table(Pactamycin_otus))
colSums(Pactamycin_df)

# Then, apply a function to count the number of samples with non-zero counts for Pactamycin
num_samples <- Pactamycin_df %>%
  dplyr::summarise(across(everything(), ~sum(.x > 0)))

# Finally, sum all the values in num_samples to get the total number of samples with non-zero counts for Pactamycin
num_samples_with_Pactamycin <- sum(num_samples)



### ANother way to do UPset plots

# Extract the OTU table and convert to a presence/absence matrix
otu_matrix <- phyloseq::otu_table(data_AMR.ps)
pa_matrix <- ifelse(otu_matrix > 0, 1, 0)

# Convert the matrix to a data frame in the wide format
df <- as.data.frame(t(pa_matrix))  # transpose the matrix before converting to data frame

# Preserve the sample names as a separate column
df$Sample <- rownames(df)

# Extract the sample metadata
sample_data <- data.frame(Sample = sample_names(data_AMR.ps), LibraryType = sample_data(data_AMR.ps)$LIbraryType)

# Add the LibraryType variable to the data frame
df <- inner_join(df, sample_data, by = "Sample")

# Convert to a long format data frame suitable for ggupset
df_upset <- df %>%
  pivot_longer(-c(Sample, LibraryType), names_to = "OTU", values_to = "Present") %>%
  mutate(Present = ifelse(Present == 1, TRUE, FALSE))

# Filter the data frame to include only the rows where Present == TRUE
df_upset <- df_upset %>% filter(Present == TRUE)

# Display the first few rows of the data frame
head(df_upset)
# Plot
upset_plot <- ggplot(df_upset, aes(x = distinct(LibraryType), fill = Present)) +
  geom_bar() +
  scale_x_upset() +
  theme(axis.text.x = element_text(angle = 90))

# Print the plot
print(upset_plot)