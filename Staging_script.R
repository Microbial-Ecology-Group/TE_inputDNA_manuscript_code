# Staging script
source("scripts/R_load_packages.R")

source("scripts/AMR_1_load_data.R")
group.ps

source("scripts/AMR_1_count_stats.R")

source("scripts/AMR_2_alpha_diversity.R")

source("scripts/AMR_3_rel_abund.R")

source("scripts/AMR_4_betadiversity.R")

source("scripts/AMR_5_upset.R")

source("scripts/AMR_6_RDA.R")

# Shotgun results
source("scripts/Shotgun_load_data.R")
shotgun_AMR_group.ps

# Combined relative abundance
source("scripts/Final_relative_abundance.R")