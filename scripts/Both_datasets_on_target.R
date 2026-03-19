library(phyloseq)
library(ggplot2)

group_palette <- c(
  "Human Wastewater" = "#7cb9cb",  # Deep blue
  "Porcine Feces" = "#FB9A99",     # Vibrant orange
  "Avian Feces" = "#33a02c",       # Green
  "Bovine Feces" = "#b15928",      # Brown
  "Canine Feces" = "#6a3d9a"       # Dark purple
)

# Extract sample data as data frames
df_te <- data.frame(sample_data(group.ps))
df_shotgun <- data.frame(sample_data(shotgun_AMR_group.ps))

# =============================================================================
# CALCULATE ON-TARGET PROPORTION
# =============================================================================

# Raw_reads * 2 because paired-end reads = 2 sequences per read pair
df_te$OnTarget <- df_te$AMR_hits / (df_te$Raw_reads * 2)
df_shotgun$OnTarget <- df_shotgun$AMR_hits / (df_shotgun$Raw_reads * 2)

# Split Part A and Part B from group.ps
df_partB <- df_te[df_te$SampleGroup == "Canine fecal", ]
df_partA <- df_te[df_te$SampleGroup != "Canine fecal", ]

# =============================================================================
# 1. SUMMARY STATISTICS
# =============================================================================

# Summary table
summary_table <- data.frame(
  Dataset = c("Part A (TE)", "Part B (TE)", "Part A (Shotgun)"),
  N = c(nrow(df_partA), nrow(df_partB), nrow(df_shotgun)),
  Mean = c(mean(df_partA$OnTarget), mean(df_partB$OnTarget), mean(df_shotgun$OnTarget)),
  SD = c(sd(df_partA$OnTarget), sd(df_partB$OnTarget), sd(df_shotgun$OnTarget)),
  Min = c(min(df_partA$OnTarget), min(df_partB$OnTarget), min(df_shotgun$OnTarget)),
  Max = c(max(df_partA$OnTarget), max(df_partB$OnTarget), max(df_shotgun$OnTarget))
)

# Format as percentages for easier interpretation
summary_table_pct <- summary_table
summary_table_pct[, 3:6] <- round(summary_table_pct[, 3:6] * 100, 4)
colnames(summary_table_pct)[3:6] <- c("Mean (%)", "SD (%)", "Min (%)", "Max (%)")
print(summary_table_pct)

# =============================================================================
# 2. STATISTICAL TESTS - Part A (TE)
# =============================================================================

cat("\n=== Part A (TE) Statistical Tests ===\n\n")

# By SampleGroup
cat("Kruskal-Wallis: OnTarget ~ SampleGroup\n")
print(kruskal.test(OnTarget ~ SampleGroup, data = df_partA))
pairwise.wilcox.test(df_partA$OnTarget ,df_partA$SampleGroup, p.adjust.method = "BH")

# By LibraryType
cat("\nKruskal-Wallis: OnTarget ~ LIbraryType\n")
print(kruskal.test(OnTarget ~ LIbraryType, data = df_partA))

# By Kit
cat("\nKruskal-Wallis: OnTarget ~ Kit\n")
print(kruskal.test(OnTarget ~ Kit, data = df_partA))

# By InputDNA
cat("\nKruskal-Wallis: OnTarget ~ InputDNA\n")
print(kruskal.test(OnTarget ~ InputDNA, data = df_partA))

# =============================================================================
# 3. STATISTICAL TESTS - Part A (Shotgun)
# =============================================================================

cat("\n=== Part A (Shotgun) Statistical Tests ===\n\n")

# By SampleGroup
cat("Kruskal-Wallis: OnTarget ~ SampleGroup\n")
print(kruskal.test(OnTarget ~ SampleGroup, data = df_shotgun))

pairwise.wilcox.test(df_shotgun$OnTarget ,df_shotgun$SampleGroup, p.adjust.method = "BH")

# =============================================================================
# 4. COMPARE Part A (TE) vs Part A (Shotgun)
# =============================================================================

cat("\n=== Part A: TE vs Shotgun Comparison ===\n\n")

# Add method label and combine
df_partA$Method <- "TE"
df_shotgun$Method <- "Shotgun"

df_compare <- rbind(
  df_partA[, c("OnTarget", "Method", "SampleGroup")],
  df_shotgun[, c("OnTarget", "Method", "SampleGroup")]
)

# Wilcoxon test
cat("Wilcoxon rank-sum test: OnTarget ~ Method (TE vs Shotgun)\n")
print(wilcox.test(OnTarget ~ Method, data = df_compare))

# =============================================================================
# 5. STATISTICAL TESTS - Part B
# =============================================================================

cat("\n=== Part B Statistical Tests ===\n\n")

cat("Kruskal-Wallis: OnTarget ~ InputDNA\n")
print(kruskal.test(OnTarget ~ InputDNA, data = df_partB))


# =============================================================================
# Prepare data with all needed columns
# =============================================================================
library(ggplot2)
library(patchwork)  # or use cowplot::plot_grid()

df_partA$Dataset <- "Part A (TE)"
df_partB$Dataset <- "Part B (TE)"
df_shotgun$Dataset <- "Part A (Shotgun)"

# For Part A TE and Shotgun combined (Panels A and B)
df_partA_all <- rbind(
  df_partA[, c("OnTarget", "Dataset", "SampleGroup", "Kit", "InputDNA")],
  df_shotgun[, c("OnTarget", "Dataset", "SampleGroup", "Kit", "InputDNA")]
)
df_partA_all$OnTarget_pct <- df_partA_all$OnTarget * 100
df_partA_all$InputDNA <- factor(df_partA_all$InputDNA)

df_partA_all$Dataset <- factor(df_partA_all$Dataset, 
                              levels = c("Part A (TE)", "Part B (TE)", "Part A (Shotgun)"))

# For all datasets (Panel A)
df_plot_all <- rbind(
  df_partA[, c("OnTarget", "Dataset", "SampleGroup", "Kit", "InputDNA")],
  df_partB[, c("OnTarget", "Dataset", "SampleGroup", "Kit", "InputDNA")],
  df_shotgun[, c("OnTarget", "Dataset", "SampleGroup", "Kit", "InputDNA")]
)
df_plot_all$OnTarget_pct <- df_plot_all$OnTarget * 100
df_plot_all$Dataset <- factor(df_plot_all$Dataset, 
                              levels = c("Part A (TE)", "Part B (TE)", "Part A (Shotgun)"))

# For Part B only (Panel C)
df_partB$OnTarget_pct <- df_partB$OnTarget * 100
df_partB$InputDNA <- factor(df_partB$InputDNA)

# =============================================================================
# Panel A: On-target by Dataset, points shaped by Kit
# =============================================================================

panelA <- ggplot(df_plot_all, aes(x = Dataset, y = OnTarget_pct, fill = Dataset)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  geom_jitter(aes(shape = Kit), width = 0.2, alpha = 0.6, size = 2) +
  scale_fill_manual(values = c("Part A (TE)" = "#1b9e77", 
                               "Part B (TE)" = "#7570b3", 
                               "Part A (Shotgun)" = "#d95f02")) +
  scale_shape_manual(values = c(16, 17, 15,12)) +  # circle, triangle, square
  labs(x = "", y = "On-target proportion (%)", shape = "Kit") +
  theme_bw() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  guides(fill = "none")


panelA

# panelA <- ggplot(df_plot_all, aes(x = Dataset, y = OnTarget_pct, fill = factor(InputDNA))) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.75)) +
#   geom_point(aes(shape = Kit), 
#              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
#              alpha = 0.6, size = 2) +
#   scale_fill_brewer(palette = "Set2", name = "Input DNA (ng)") +
#   scale_shape_manual(values = c(16, 17, 15, 12), name = "Kit") +
#   labs(x = "", y = "On-target proportion (%)") +
#   theme_bw() +
#   theme(legend.position = "right",
#         axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
# 
# panelA


# =============================================================================
# Panel B: Part A (TE + Shotgun) by SampleGroup, dodged by Dataset, shaped by InputDNA
# =============================================================================

panelB <- ggplot(df_partA_all, aes(x = SampleGroup, y = OnTarget_pct, fill = Dataset)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.75)) +
  geom_point(aes(shape = Kit, group = Dataset),
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
             alpha = 0.6, size = 2) +
  scale_fill_manual(values = c("Part A (TE)" = "#1b9e77", 
                               "Part A (Shotgun)" = "#d95f02")) +
  scale_shape_manual(values = c(16, 17, 15, 12), name = "Kit") +
  labs(x = "Sample Group", y = "On-target proportion (%)", 
       fill = "Dataset") +
  theme_bw() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

panelB
# =============================================================================
# Panel C: Part B only, by InputDNA
# =============================================================================

panelC <- ggplot(df_partB, aes(x = InputDNA, y = OnTarget_pct, fill = InputDNA)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, varwidth = TRUE) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 2) +
  scale_fill_manual(values = c("#7570b3", "#beaed4")) +  # shades of purple for Part B
  labs(x = "Input DNA (ng)", y = "On-target proportion (%)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10))

panelC

# Combine Part A (TE only) and Part B
df_partAB_TE <- rbind(
  df_partA[, c("OnTarget", "Dataset", "SampleGroup", "Kit", "InputDNA")],
  df_partB[, c("OnTarget", "Dataset", "SampleGroup", "Kit", "InputDNA")]
)
df_partAB_TE$OnTarget_pct <- df_partAB_TE$OnTarget * 100
df_partAB_TE$InputDNA <- factor(df_partAB_TE$InputDNA)

panelC <- ggplot(df_partAB_TE, aes(x = SampleGroup, y = OnTarget_pct, fill = InputDNA)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.75)) +
  geom_point(aes(shape = Kit, group = InputDNA),
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
             alpha = 0.6, size = 2) +
  scale_fill_brewer(palette = "Set2", name = "Input DNA (ng)") +
  scale_shape_manual(values = c(16, 17, 15, 12), name = "Kit") +
  labs(x = "Sample Group", y = "On-target proportion (%)") +
  theme_bw() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

panelC


# =============================================================================
# Combine panels
# =============================================================================

# Option 1: Using patchwork
combined <- panelA + panelB + panelC +
  plot_layout(ncol = 2, widths = c(1, 1, 2)) +
  plot_annotation(tag_levels = 'A')

print(combined)

# Save
ggsave("ontarget_3panel.png", combined, width = 14, height = 5, dpi = 300)

# =============================================================================
# Option 2: Stacked layout (if horizontal is too cramped)
# =============================================================================

combined_stacked <- panelA / (panelB | panelC) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = 'A')

print(combined_stacked)

ggsave("ontarget_3panel_stacked.png", combined_stacked, width = 12, height = 8, dpi = 300)
# =============================================================================
# 7. FOLD-CHANGE CALCULATION (TE vs Shotgun)
# =============================================================================

cat("\n=== Fold-change: TE vs Shotgun ===\n")
fold_change <- mean(df_partA$OnTarget) / mean(df_shotgun$OnTarget)
cat("Mean on-target (TE):", round(mean(df_partA$OnTarget) * 100, 4), "%\n")
cat("Mean on-target (Shotgun):", round(mean(df_shotgun$OnTarget) * 100, 4), "%\n")
cat("Fold enrichment (TE/Shotgun):", round(fold_change, 1), "x\n")

# Supplemental: On-target by InputDNA (Part A only)
p_input <- ggplot(df_partA, aes(x = factor(InputDNA), y = OnTarget, fill = factor(InputDNA))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(x = "DNA Input Amount (ng)", y = "On-target proportion (%)",
       title = "On-target proportion by DNA input (Part A)") +
  theme_bw() +
  theme(legend.position = "none")

print(p_input)
              