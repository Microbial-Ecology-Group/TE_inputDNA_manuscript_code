# Count stats
colnames(sample_data(group.ps))
colnames(sample_data(shotgun_AMR_group.ps))


# Extract sample data as data frames
df_te <- data.frame(sample_data(group.ps))
df_shotgun <- data.frame(sample_data(shotgun_AMR_group.ps))

# Split Part A and Part B from group.ps
df_partB <- df_te[df_te$SampleGroup == "Canine fecal", ]
df_partA <- df_te[df_te$SampleGroup != "Canine fecal", ]

# =============================================================================
# 1. SUMMARY STATISTICS
# =============================================================================

# Part A (TE)
cat("=== Part A (TE) ===\n")
cat("N:", nrow(df_partA), "\n")
cat("Mean Raw_reads:", mean(df_partA$Raw_reads), "\n")
cat("Range:", min(df_partA$Raw_reads), "-", max(df_partA$Raw_reads), "\n\n")

# Part B (TE)
cat("=== Part B (TE - Canine) ===\n")
cat("N:", nrow(df_partB), "\n")
cat("Mean Raw_reads:", mean(df_partB$Raw_reads), "\n")
cat("Range:", min(df_partB$Raw_reads), "-", max(df_partB$Raw_reads), "\n\n")

# Part A (Shotgun)
cat("=== Part A (Shotgun) ===\n")
cat("N:", nrow(df_shotgun), "\n")
cat("Mean Raw_reads:", mean(df_shotgun$Raw_reads), "\n")
cat("Range:", min(df_shotgun$Raw_reads), "-", max(df_shotgun$Raw_reads), "\n\n")

# Nice summary table
summary_table <- data.frame(
  Dataset = c("Part A (TE)", "Part B (TE)", "Part A (Shotgun)"),
  N = c(nrow(df_partA), nrow(df_partB), nrow(df_shotgun)),
  Mean = c(mean(df_partA$Raw_reads), mean(df_partB$Raw_reads), mean(df_shotgun$Raw_reads)),
  Min = c(min(df_partA$Raw_reads), min(df_partB$Raw_reads), min(df_shotgun$Raw_reads)),
  Max = c(max(df_partA$Raw_reads), max(df_partB$Raw_reads), max(df_shotgun$Raw_reads))
)
print(summary_table)

# =============================================================================
# 2. STATISTICAL TESTS - Part A (TE)
# =============================================================================

cat("\n=== Part A (TE) Statistical Tests ===\n\n")

# By SampleGroup
cat("Kruskal-Wallis: Raw_reads ~ SampleGroup\n")
kruskal.test(Raw_reads ~ SampleGroup, data = df_partA)

# By LibraryType
cat("\nKruskal-Wallis: Raw_reads ~ LIbraryType\n")
kruskal.test(Raw_reads ~ LIbraryType, data = df_partA)

# By Kit
cat("\nKruskal-Wallis: Raw_reads ~ Kit\n")
kruskal.test(Raw_reads ~ Kit, data = df_partA)

# By InputDNA
cat("\nKruskal-Wallis: Raw_reads ~ InputDNA\n")
kruskal.test(Raw_reads ~ InputDNA, data = df_partA)

# =============================================================================
# 3. STATISTICAL TESTS - Part A (Shotgun)
# =============================================================================

cat("\n=== Part A (Shotgun) Statistical Tests ===\n\n")

# By SampleGroup
cat("Kruskal-Wallis: Raw_reads ~ SampleGroup\n")
kruskal.test(Raw_reads ~ SampleGroup, data = df_shotgun)

# =============================================================================
# 4. COMPARE Part A (TE) vs Part A (Shotgun)
# =============================================================================

cat("\n=== Part A: TE vs Shotgun Comparison ===\n\n")

# Add method label and combine
df_partA$Method <- "TE"
df_shotgun$Method <- "Shotgun"

# Get common columns for comparison
df_compare <- rbind(
  df_partA[, c("Raw_reads", "Method")],
  df_shotgun[, c("Raw_reads", "Method")]
)

# Wilcoxon test
cat("Wilcoxon rank-sum test: Raw_reads ~ Method (TE vs Shotgun)\n")
wilcox.test(Raw_reads ~ Method, data = df_compare)

# =============================================================================
# 5. STATISTICAL TESTS - Part B
# =============================================================================

cat("\n=== Part B Statistical Tests ===\n\n")

cat("\nKruskal-Wallis: Raw_reads ~ InputDNA\n")
kruskal.test(Raw_reads ~ InputDNA, data = df_partB)

