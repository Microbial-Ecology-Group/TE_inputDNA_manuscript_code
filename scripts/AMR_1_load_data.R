#### TE DATA ####
# read in AMR count matrix
amr_data <- read.table('MEG_v3_results/deduped_SNPconfirmed_AMR_analytic_matrix_AgilentInput.csv', header=T, row.names=1, sep=',', quote = "")
# convert this into an 'otu_table' format required for phyloseq
amr_data <- otu_table(amr_data, taxa_are_rows = T)

#read in gene annotations
annotations <- read.table('MEG_v3_results/megares_annotations_v3.00.csv', header=T, row.names=1, sep=",", quote = "")
#convert this into a 'taxonomy table' format required for phyloseq
annotations <- phyloseq::tax_table(as.matrix(annotations))

# read in TE metadata
amr_metadata <- read.table('Agilent_input_metadata.csv', header=T, sep=',', row.names = 1, quote = "")
# convert to 'sample_data' format for phyloseq
amr_metadata <- sample_data(amr_metadata)

# merge the annotations, the count matrix, and metadata into a phyloseq object
amr.ps <- merge_phyloseq(amr_data, annotations, amr_metadata)

# Loaded object
amr.ps # 110 samples

### need to split up the genes needing SNP confirmation:
data_SNPconfirm <- subset_taxa(amr.ps, snp=="RequiresSNPConfirmation", T) # 395 of the 4037 genes
data_noSNP <- subset_taxa(amr.ps, snp!="RequiresSNPConfirmation")
data_noSNP # 3642 of the 4037 genes

# We keep the SNP confirmed counts anyway
data_AMR.ps <- amr.ps

# # some QC checks
sum(sample_sums(data_noSNP)) # 89298350, 93103894
min(sample_sums(data_noSNP)) # 123310
max(sample_sums(data_noSNP)) # 1336239
sort(sample_sums(data_noSNP)) # 

sum(sample_sums(data_SNPconfirm)) # 4230878
min(sample_sums(data_SNPconfirm)) # 4157
max(sample_sums(data_SNPconfirm)) # 77711
sort(sample_sums(data_SNPconfirm)) 




# Look at sorted sample sumes
sort(sample_sums(data_AMR.ps)) # 539, 552, 3641, 5177, 5294, 6591, 16,532, 16,886...

sort(taxa_sums(data_AMR.ps)) # 539, 552, 3641, 5177, 5294, 6591, 16,532, 16,886...


# # for now, going to cut off those below 3921 but I might be inclined to do 16,886
# 
# data_noSNP <- prune_samples(sample_sums(data_noSNP) > 3500, data_noSNP)
# data_noSNP # lost two samples and now have 2618 taxa
# any(taxa_sums(data_noSNP)==0) # no, good



## Make new data columns

sample_data(data_AMR.ps)$AMR_hits <- sample_sums(data_AMR.ps)



# Check SNP removed samples
data_AMR.ps # 111 samples



### Make metadata file to use later ####
AMR_mapfile <- data.frame(sample_data(data_AMR.ps))
AMR_mapfile$Dilution <- as.factor(AMR_mapfile$Dilution)
AMR_mapfile$InputDNA <- as.factor(AMR_mapfile$InputDNA)
AMR_mapfile$LIbraryType <- as.factor(AMR_mapfile$LIbraryType)
AMR_mapfile$SampleGroup <- as.factor(AMR_mapfile$SampleGroup)
AMR_mapfile$OriginalSample <- as.factor(AMR_mapfile$OriginalSample)

AMR_mapfile$Scaled_Raw_reads <- scale(AMR_mapfile$Raw_reads)

AMR_mapfile$Dilution <- factor(AMR_mapfile$Dilution, levels = c("Standard","1:2 dilution",
                                                                      "1:3 dilution", "1:4 dilution","1:5 dilution"))

# Check which variables have any missing data
sapply(AMR_mapfile, function(x) any(is.na(x)))

group.ps <- tax_glom(data_AMR.ps, taxrank = "group")

