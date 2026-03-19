# Manuscript figures
library(cowplot)
library(UpSetR)
library(tidyverse)
library(MicrobiotaProcess)

## Figure 1 - Study design (Biorender) ####

## Figure 2  ####
# Consider making a panel with AMR hit by dilution, then smaller boxplots comparing by species
boxplot(AMR_mapfile$AMR_hits ~ AMR_mapfile$Dilution)



## Figure 2 #####
richness_plot <- ggplot(alpha_div_meta[alpha_div_meta$SampleGroup != "Canine fecal", ], aes(x= SampleGroup, y= Observed, fill = LIbraryType, colour = LIbraryType)) +
  theme_bw() + 
  labs(y= "Observed ARGs") +
  geom_boxplot(alpha=0.4) +
  geom_point(position = position_dodge(width = 0.8)) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size =24, colour = "white"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 28),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 28),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

shannon_plot <- ggplot(alpha_div_meta[alpha_div_meta$SampleGroup != "Canine fecal", ], aes(x= SampleGroup, y= Shannon, fill = LIbraryType, colour = LIbraryType)) +
  theme_bw() + 
  labs(y= "Shannon's index") +
  geom_boxplot(alpha=0.4) +
  geom_point(position = position_dodge(width = 0.8)) +
  theme(legend.position = "right",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "black", size = 1.0),
        strip.text = element_text(size =24, colour = "white"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 28),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 28),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  shannon_plot + theme(legend.box.margin = margin(0, 0, 0, 12), legend.text = element_text(size = 14),
                       legend.title = element_text(size = 16))
)

# Combine the plots
combined_plot_row <- plot_grid(richness_plot, shannon_plot + theme(legend.position="none"),  
                               align = 'vh',
                               hjust = -1,
                               nrow = 1,
                           labels = "AUTO")


# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
plot_grid(combined_plot_row, legend, rel_widths = c(3, .4))


## Upset Figure ####

# Kit
DatasetA_group.ps <- subset_samples(ps_group, SampleGroup != "Canine fecal")
group_libkit_upsetda_edit <- get_upset(DatasetA_group.ps, factorNames="LIbraryType") ## ASV


group_level_upset <- upset(group_libkit_upsetda_edit, sets=c("XT","XT_HS"),
                          sets.bar.color = c("red", "blue"),text.scale = 2,
                          order.by = "freq", empty.intersections = "on")

group_level_upset


# Kit
# DatasetA_class.ps <- subset_samples(ps_class, SampleGroup != "Canine fecal")
# class_libkit_upsetda_edit <- get_upset(DatasetA_class.ps, factorNames="LIbraryType") ## ASV
# 
# 
# class_level_upset <- upset(class_libkit_upsetda_edit, sets=c("XT","XT_HS"),
#                            sets.bar.color = c("red", "blue"),text.scale = 2,
#                            order.by = "freq", empty.intersections = "on")
# 
# class_level_upset

## Figure 3 - Diversity #####


### Figure 4 - Relative abundance dendrogram ####


### Figure ANCOMBC ####
fig_ANCOMBC_grid


## Supplemental ####
### Sup 1 - Ordination all samples, by species and kit ####

AllSamples.plt_ord_by_Group_library <- plot_ordination(data_AMR.ps.css, AllSamples.data_AMR.ps.css.ord, color = "SampleGroup", shape = "LIbraryType") +
  theme_bw() +
  labs(title ="Sample type") +
  stat_ellipse(aes(fill= SampleGroup), geom="polygon", alpha = 0.25) +
  theme(legend.position = "right",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_rect(fill= "grey91", linewidth = 1.0),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.7),
        plot.title = element_text(size = 30),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
AllSamples.plt_ord_by_Group_library
