library(dplyr)
library(purrr)

# Function to filter groups and run the pairwise Wilcoxon test
run_tests_by_group <- function(data, test_column, group_column, split_column) {
  
  # Step 1: Filter out groups that don't have at least 2 levels in group_column
  valid_data <- data %>%
    group_by_at(split_column) %>%
    filter(n_distinct(.data[[group_column]]) > 1) %>%
    ungroup()
  
  # Step 2: Split the valid data by split_column and apply the pairwise Wilcoxon test
  results <- valid_data %>%
    group_by_at(split_column) %>%
    group_split() %>%
    map(~ pairwise.wilcox.test(
      .x[[test_column]],   # Column to test (e.g., Observed richness)
      .x[[group_column]],  # Column to group by (e.g., LIbraryType)
      p.adjust.method = "BH" # Adjust p-values with Benjamini-Hochberg
    ))
  
  # Step 3: Assign names to the list of results based on split_column values
  names(results) <- unique(valid_data[[split_column]])
  
  return(results)
}

# Example usage
valid_results <- run_tests_by_group(
  alpha_div_meta,          # Data frame
  test_column = "Observed", # Column to test (e.g., Observed richness)
  group_column = "LIbraryType", # Column to group by (e.g., LIbraryType)
  split_column = "SampleGroup"  # Column to split the data by (e.g., SampleGroup)
)

# Access the results for a specific SampleGroup
valid_results
