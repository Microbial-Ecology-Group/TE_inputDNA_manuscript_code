# Code for comparing raw reads and AMR hits between TE and shotgun
library(dplyr)
sample_hits.df <- read.table("MEG_v3_results/Agilent_TE_vs_shotgun.csv", header = TRUE, sep = ",")

sample_hits.df$InputDNA <- factor(sample_hits.df$InputDNA, levels = c("Shotgun", "50ng" ,"200ng", "400ng", "800ng"))

sample_hits.df <- sample_hits.df %>% 
  mutate(    Percent_on_target = round( (AMR_hits / (Raw_reads * 2)) * 100, 2 )   )

## Comparing Shotgun to TE ####
group_palette <- c( "navyblue","gold2")
kit_palette <- c("navyblue", "darkgreen", "maroon4")

# Raw reads by SampleGroup
raw_reads_plot <- plot_nested_wilcox_letters(
  sample_hits.df %>% filter(SampleGroup != "Canine fecal" ), #
  yvar = "Raw_reads",
  ylab = "Raw paired reads",
  palette = group_palette,
  test_var = "Type",
  group_var = "SampleGroup",
  compare_scope = "global"
)

png("figures/Fig_Raw_reads_Shotgun_vs_TE_by_sample_group.png", width = 1450, height = 1125)
raw_reads_plot
dev.off()


# AMR hits by SampleGroup
AMR_hits_plot <- plot_nested_wilcox_letters(
  sample_hits.df %>% filter(SampleGroup != "Canine fecal" ),
  yvar = "AMR_hits",
  ylab = "AMR_hits",
  palette = group_palette,
  test_var = "Type",
  group_var = "SampleGroup",
  compare_scope = "global"
)

# Percent_on_target by SampleGroup
Percent_on_target_plot <- plot_nested_wilcox_letters(
  sample_hits.df %>% filter(SampleGroup != "Canine fecal" ),
  yvar = "Percent_on_target",
  ylab = "Percent on target",
  palette = group_palette,
  test_var = "Type",
  group_var = "SampleGroup",
  compare_scope = "global"
)

# Plots comparing Type
Observed_plot <- plot_nested_wilcox_letters(
  sample_hits.df %>% filter(SampleGroup != "Canine fecal" ),
  yvar = "Observed",
  ylab = "Observed",
  palette = group_palette,
  test_var = "Type",
  group_var = "SampleGroup",
  compare_scope = "global"
)

Shannon_plot <- plot_nested_wilcox_letters(
  sample_hits.df %>% filter(SampleGroup != "Canine fecal" ),
  yvar = "Shannon",
  ylab = "Shannon",
  palette = group_palette,
  test_var = "Type",
  group_var = "SampleGroup",
  compare_scope = "global"
)

# Plots comparing Kit
Observed_plot <- plot_nested_wilcox_letters(
  sample_hits.df %>% filter(SampleGroup != "Canine fecal"  ),
  yvar = "Observed",
  ylab = "Observed",
  palette = kit_palette,
  test_var = "Kit",
  group_var = "SampleGroup",
  compare_scope = "within"
)

Shannon_plot <- plot_nested_wilcox_letters(
  sample_hits.df %>% filter(SampleGroup != "Canine fecal" ),
  yvar = "Shannon",
  ylab = "Shannon",
  palette = kit_palette,
  test_var = "Kit",
  group_var = "SampleGroup",
  compare_scope = "within"
)


diversity_grid_by_kit <- plot_grid(
  Observed_plot + theme(legend.position = "none"),
  Shannon_plot+ theme(legend.position = "none"),
  cowplot::get_legend(Shannon_plot),
  ncol = 3,
  labels = c("A", "B"),
  rel_widths = c(1,1,.5)
)


# Save the final plot
#ggsave("../Writing/Final_figures/Fig3_diversity_boxplot_PartA_by_kit.jpg", plot = diversity_grid_by_kit, width = 20, height = 11, dpi = 300)

png("../Writing/Final_figures/Fig3_diversity_boxplot_PartA_by_kit.png", width = 2400, height = 1400, res = 100)
diversity_grid_by_kit
dev.off()

#
##
# Comparing InputDNA by sample group ####
##
#

inputDNA_palette <- c(
  "navyblue",      # keeps your original purple
  #"#FCE620",      # vivid yellow
  "gold2",          # classic gold
  "#F7981C",      # rich orange
  "brown4"
)

named_inputDNA_palette <- c(
  "shotgun" = "navyblue",      # keeps your original purple
  "50ng" = "#FCE620",      # vivid yellow
  "200ng" = "gold2",          # classic gold
  "400ng" = "#F7981C",      # rich orange
  "800ng" = "brown4"
)

## Avian fecal Raw_reads
Avian_raw_read_plot <- plot_wilcox_letters(
  sample_hits.df %>% filter(SampleGroup == "Avian fecal"),
  yvar = "Raw_reads",
  ylab = "Raw paired reads",
  palette = inputDNA_palette,
  test_var = "InputDNA"
)
 
## Canine fecal Raw_reads
Canine_raw_read_plot <- plot_wilcox_letters(
  sample_hits.df %>% filter(SampleGroup == "Canine fecal"),
  yvar = "Raw_reads",
  ylab = "Raw paired reads",
  palette = named_inputDNA_palette,
  test_var = "InputDNA"
)

## Bovine fecal Raw_reads
Bovine_raw_read_plot <- plot_wilcox_letters(
  sample_hits.df %>% filter(SampleGroup == "Bovine fecal"),
  yvar = "Raw_reads",
  ylab = "Raw paired reads",
  palette = inputDNA_palette,
  test_var = "InputDNA"
)

## Porcine fecal Raw_reads
Porcine_raw_read_plot <- plot_wilcox_letters(
  sample_hits.df %>% filter(SampleGroup == "Porcine fecal"),
  yvar = "Raw_reads",
  ylab = "Raw paired reads",
  palette = inputDNA_palette,
  test_var = "InputDNA"
)

## Human wastewater Raw_reads
Human_wastewater_raw_read_plot <- plot_wilcox_letters(
  sample_hits.df %>% filter(SampleGroup == "Human wastewater"),
  yvar = "Raw_reads",
  ylab = "Raw paired reads",
  palette = inputDNA_palette,
  test_var = "InputDNA"
)

## Avian fecal Percent_on_target
Avian_Percent_on_target_plot <- plot_wilcox_letters(
  sample_hits.df %>% filter(SampleGroup == "Avian fecal"),
  yvar = "Percent_on_target",
  ylab = "Percent on-target",
  palette = inputDNA_palette,
  test_var = "InputDNA"
)

## Canine fecal Percent_on_target
Canine_Percent_on_target_plot <- plot_wilcox_letters(
  sample_hits.df %>% filter(SampleGroup == "Canine fecal"),
  yvar = "Percent_on_target",
  ylab = "Percent on-target",
  palette = named_inputDNA_palette,
  test_var = "InputDNA"
)

## Bovine fecal Percent_on_target
Bovine_Percent_on_target_plot <- plot_wilcox_letters(
  sample_hits.df %>% filter(SampleGroup == "Bovine fecal"),
  yvar = "Percent_on_target",
  ylab = "Percent on-target",
  palette = inputDNA_palette,
  test_var = "InputDNA"
)

## Porcine fecal Percent_on_target
Porcine_Percent_on_target_plot <- plot_wilcox_letters(
  sample_hits.df %>% filter(SampleGroup == "Porcine fecal"),
  yvar = "Percent_on_target",
  ylab = "Percent on-target",
  palette = inputDNA_palette,
  test_var = "InputDNA"
)

## Human wastewater Percent_on_target
Human_wastewater_Percent_on_target_plot <- plot_wilcox_letters(
  sample_hits.df %>% filter(SampleGroup == "Human wastewater"),
  yvar = "Percent_on_target",
  ylab = "Percent on-target",
  palette = inputDNA_palette,
  test_var = "InputDNA"
)

## Make cowplot of count comparison by input DNA
library(cowplot)
library(ggplot2)

## ──────────────────────────────────────────────────────────────
## 0 •    names here must match your actual object names
##        (change if needed)
## ──────────────────────────────────────────────────────────────
raw_list <- list(
  Avian_raw_read_plot,
  Bovine_raw_read_plot,
  Porcine_raw_read_plot,
  Human_wastewater_raw_read_plot
)

amr_list <- list(
  Avian_Percent_on_target_plot,
  Bovine_Percent_on_target_plot,
  Porcine_Percent_on_target_plot,
  Human_wastewater_Percent_on_target_plot
)

final_plot <- plot_grid(  Avian_raw_read_plot,
                          Bovine_raw_read_plot,
                          Porcine_raw_read_plot,
                          Human_wastewater_raw_read_plot,
                          Avian_Percent_on_target_plot,
                          Bovine_Percent_on_target_plot,
                          Porcine_Percent_on_target_plot,
                          Human_wastewater_Percent_on_target_plot,
                        nrow = 2, align = "hv")

canine_plot <- plot_grid( Canine_raw_read_plot,Canine_Percent_on_target_plot,nrow=2, align = "hv")

# (raw_list and amr_list already created; each has 4 ggplots)

# all_panels <- c(raw_list, amr_list)           # eight plots in order
# 
# shared_legend <- get_legend(raw_list[[1]] +
#                               theme(legend.position = "right",
#                                     legend.title = element_text(size = 30),
#                                     legend.text = element_text(size = 24)))
# 
# # remove legends and keep default margins
# all_panels <- lapply(all_panels, function(p)
#   p + theme(legend.position = "none",
#             plot.margin = margin(5.5,5.5,5.5,5.5)))
# 
# full_grid <- plot_grid(plotlist = all_panels,
#                        nrow   = 2,
#                        align  = "hv",
#                        labels = LETTERS[1:8],
#                        label_size = 18)
# 
# final_plot <- plot_grid(full_grid, shared_legend,
#                         ncol = 2, rel_widths = c(1, .08))

#ggsave("Raw_vs_AMR_2x4_grid.pdf", final_plot, width = 14, height = 7)

get_legend <- function(p)
  cowplot::get_legend(
    p + theme(
      legend.position = "right",
      legend.title    = element_text(size = 30),
      legend.text     = element_text(size = 24)
    )
  )

shared_legend <- get_legend(raw_list[[1]])

# Canine plot
Canine_legend <- get_legend(Canine_raw_read_plot) 

############################################################################
##  2 • strip legends, keep identical margins
############################################################################
fix_panel <- function(p)
  p + theme(
    legend.position = "none",
    plot.margin     = margin(5.5, 5.5, 5.5, 5.5)  # identical on all panels
  )

all_panels <- c(raw_list, amr_list)
all_panels <- lapply(all_panels, fix_panel)


## Canine plot
canine_panels <- c(Canine_raw_read_plot,Canine_Percent_on_target_plot)
canine_panels <- lapply(canine_panels, fix_panel)

############################################################################
##  3 • build the 2 × 4 grid (A–H)
############################################################################
full_grid <- cowplot::plot_grid(
  plotlist   = all_panels,
  nrow       = 2,
  align      = "hv",
  axis       = "tblr",
  labels     = LETTERS[1:8],
  label_size = 30,
  label_fontface = "bold"
)

canine_grid <- cowplot::plot_grid(
  plotlist   = canine_panels,
  nrow       = 2,
  align      = "hv",
  axis       = "tblr",
  labels     = LETTERS[1:8],
  label_size = 30,
  label_fontface = "bold"
)

############################################################################
##  4 • add an external header strip with column names
############################################################################
col_titles <- c("Avian fecal", "Bovine fecal",
                "Porcine fecal", "Human wastewater")

header <- cowplot::ggdraw() +
  cowplot::theme_nothing() +
  purrr::imap(col_titles, \(lab, i) {
    cowplot::draw_label(
      lab,
      x        = (i - .5) / length(col_titles) + .01,  # centred over each column
      y        = .5,
      fontface = "bold",
      size     = 30,
      hjust    = .5
    )
  })

with_header <- cowplot::plot_grid(
  header, full_grid,
  ncol        = 1,
  rel_heights = c(.06, 1)           # header strip height
)


## Canine plot

canine_col_titles <- c("Canine fecal")

canine_header <- cowplot::ggdraw() +
  cowplot::theme_nothing() +
  purrr::imap(canine_col_titles, \(lab, i) {
    cowplot::draw_label(
      lab,
      x        = (i - .5) / length(canine_col_titles) + .01,  # centred over each column
      y        = .5,
      fontface = "bold",
      size     = 30,
      hjust    = .5
    )
  })

canine_with_header <- cowplot::plot_grid(
  canine_header, canine_grid,
  ncol        = 1,
  rel_heights = c(.06, 1)           # header strip height
)


############################################################################
##  5 • place legend on the right
############################################################################
final_plot <- cowplot::plot_grid(
  with_header, shared_legend,
  ncol        = 2,
  rel_widths  = c(1, .08)
)

canine_plot <- cowplot::plot_grid(
  canine_with_header, Canine_legend,
  ncol        = 2,
  rel_widths  = c(.5, .08)
)


############################################################################
##  6 • view / save
############################################################################
print(final_plot)


png("figures/Fig_Shotgun_vs_TE_facet_sample_group_byInputDNA.png", width = 2450, height = 1125)
final_plot
dev.off()


png("figures/Fig_Canine_byInputDNA.png", width = 1000, height = 1125)
canine_plot
dev.off()
