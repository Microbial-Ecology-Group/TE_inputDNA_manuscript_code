# Upset plot
library(UpSetR)
library(MicrobiotaProcess)
shotgun_AMR_group.ps
group.ps

all_AMR_group.ps <- merge_phyloseq(shotgun_AMR_group.ps,group.ps)
all_AMR_group.ps <- prune_taxa(taxa_sums(all_AMR_group.ps) > 0, all_AMR_group.ps)

# Kit
DatasetA_group.ps <- subset_samples(
  all_AMR_group.ps,
  SampleGroup != "Canine fecal"
)

DatasetA_group.ps <- subset_samples(
  all_AMR_group.ps,
  SampleGroup != "Canine fecal" )

sample_data(DatasetA_group.ps)$Type


## By InputDNA
group_InputDNA_upsetda_edit <- get_upset(DatasetA_group.ps, factorNames="InputDNA") ## ASV


group_InputDNA_upset <- upset(group_InputDNA_upsetda_edit, sets=c("Shotgun", "200ng", "400ng", "800ng"),
                           sets.bar.color = c(
                             "navyblue",      # keeps your original purple
                             "#FCE620",      # vivid yellow
                             "#F7981C",      # rich orange
                             "gold2"          # classic gold
                           ),text.scale = 2,
                           order.by = "freq")

group_InputDNA_upset

## By Type
group_Type_upsetda_edit <- get_upset(DatasetA_group.ps, factorNames="Type") ## ASV


group_Type_upset <- upset(group_Type_upsetda_edit, sets=c("Shotgun", "TE"),
                              sets.bar.color = c(
                                "gold2",      # keeps your original purple
                                "navyblue"          # classic gold
                              ),text.scale = 2,
                              order.by = "freq")

group_Type_upset
group_level_upset


