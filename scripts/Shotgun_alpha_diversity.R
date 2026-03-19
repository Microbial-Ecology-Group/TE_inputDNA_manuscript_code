shotgun_alpha_div <- estimate_richness(shotgun_AMR_group.ps, measures = c("Observed","Shannon","Simpson","InvSimpson"))
shotgun_alpha_div
shotgun_alpha_div.df <- as(sample_data(shotgun_AMR_group.ps), "data.frame")
shotgun_alpha_div_meta <- cbind(shotgun_alpha_div, shotgun_alpha_div.df)


write.csv(shotgun_alpha_div_meta, "agilent_shotgun_alpha_div_measures.csv")
