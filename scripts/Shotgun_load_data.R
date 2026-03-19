#### TE DATA ####
# read in AMR count matrix
shotgun_data <- read.table('Shotgun_results/Results/SNPconfirmed_dedup_AMR_analytic_matrix.csv', header=T, row.names=1, sep=',', quote = "")
# convert this into an 'otu_table' format required for phyloseq
shotgun_data <- otu_table(shotgun_data, taxa_are_rows = T)

#read in gene annotations
annotations <- read.table('MEG_v3_results/megares_annotations_v3.00.csv', header=T, row.names=1, sep=",", quote = "")
#convert this into a 'taxonomy table' format required for phyloseq
annotations <- phyloseq::tax_table(as.matrix(annotations))

# read in TE metadata
shotgun_metadata <- read.table('Shotgun_results/metadata_shotgun.csv', header=T, sep=',', row.names = 1, quote = "")
# convert to 'sample_data' format for phyloseq
shotgun_metadata <- sample_data(shotgun_metadata)

# merge the annotations, the count matrix, and metadata into a phyloseq object
shotgun_AMR.ps <- merge_phyloseq(shotgun_data, annotations, shotgun_metadata)

shotgun_AMR_group.ps <- tax_glom(shotgun_AMR.ps, taxrank = "group")
