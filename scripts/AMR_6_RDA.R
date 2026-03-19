# RDA

#DatasetA_group.ps <- subset_samples(ps_group, SampleGroup != "Canine fecal")
#DatasetA_group.ps <- data_AMR.ps


cleanset <- data.frame(sample_data(group.ps))
cleanset <- cleanset[, c("InputDNA","SampleGroup","Kit","Dilution","OriginalSample")]

#cleanset <- cleanset[, c("InputDNA","SampleGroup","Kit","Dilution")]


############################################### RDA MODEL ###################################################

#### CSS transformation of counts
group.ps.css <- phyloseq_transform_css(group.ps, log = F)

norm.mt<-as(otu_table(group.ps.css),"matrix")
### This matrix has rows as otus and columns as samples (all our work before used columns as otus)
## we first remove all rows (otus) that have count zero by summing per row and keeping those rows with
## sum greater than zero
norm.rs = apply(norm.mt,1,sum)
norm.mt = norm.mt[norm.rs>0,]
## We then transform the normalized matrix so that rows are samples and columns are otus (for vegan to work correctly)
norm.mt = t(norm.mt)
#norm.mt<-norm.mt[cleanset_samples,]
#Using function decostand() in vegan to transform using hellinger standardizations
hell.norm.mt <- decostand(norm.mt,"hell")
### We define the full model that includes all possible variables in the metadata file
mod1 <- rda(hell.norm.mt ~ ., data =  cleanset)
anova(mod1, by = "term", perm = 1000)

### We define the null model that only includes an intercept (1)
mod0 <- rda(hell.norm.mt ~ 1, data =  cleanset)
### We use the step function to find the best model by removing one term at a time.
### We use permutation test (perMANOVA) to do so
mod <- step(mod0, scope = formula(mod1), test = "perm", direction="both")
anova(mod, by = "term", perm = 1000)
anova(mod,perm = 1000)

apriori_mod <- rda(hell.norm.mt ~ SampleGroup + Kit + Dilution + InputDNA, data =  cleanset)
anova(apriori_mod, by = "term", perm = 1000)
anova(apriori_mod,perm = 1000)

# With OriginalSample
# Inertia Proportion Rank
# Total         0.255075   1.000000     
# Constrained   0.247247   0.969312   38
# Unconstrained 0.007828   0.030688   71

# Df Variance      F Pr(>F)    
# Model    38 0.247247 59.016  0.001 ***
#   Residual 71 0.007828   

# Model: rda(formula = hell.norm.mt ~ OriginalSample + Kit + Dilution + InputDNA, data = cleanset)
# Df Variance       F Pr(>F)    
# OriginalSample 31 0.242628 70.9910  0.001 ***
#   Kit     1 0.003170 28.7537  0.001 ***
#   Dilution        4 0.001010  2.2896  0.010 ** 
#   InputDNA        2 0.000439  1.9911  0.061 .  
# Residual       71 0.007828    




# View the anova results
mod

# Visualize
#### 2D RDA plot (this doesn't work for db-RDA, plots are meaningless there)
fig <-ordiplot(mod,type="n")
colvec <- c("black","grey51","blue","red","green") 
#pchvec <-c(10,2)
points(mod, "sites", pch=21, col=as.factor(scaleddf$var), bg=as.factor(scaleddf$var))
groupz <- sort(unique(scaleddf$var))
for(i in seq(groupz)) {ordispider(mod,  scaleddf$catvar,font=2, cex=1.5, col=i, show.groups=groupz[i])}
legv <- sort(unique(scaleddf$catvar))
#legend("right", legend = levels(scaled_metadata$Group),col=legv, title = "Sample groups",bty = "n", pt.bg = legv)
text(fig, "biplot", col="red", cex=0.9)
text(fig, "centroids", col=c("black","red","green"), cex=2, pos = 4)


# Apriori model
# We can define our own model here
mod_apriori <- rda(hell.norm.mt ~ InputDNA + SampleGroup + Kit + Dilution + OriginalSample , data =  cleanset)
anova(mod_apriori,perm = 10000)
anova(mod_apriori, by = "term", perm = 10000)

# All sample

# Df Variance      F Pr(>F)    
# Model    38 0.247247 59.016  0.001 ***
#   Residual 71 0.007828                  

# Inertia Proportion Rank
# Total         0.255075   1.000000     
# Constrained   0.247247   0.969312   38
# Unconstrained 0.007828   0.030688   71



# Df Variance        F Pr(>F)    
# InputDNA        3 0.022887  69.1971  0.001 ***
#   SampleGroup     4 0.170101 385.7195  0.001 ***
#   Kit     1 0.003187  28.9038  0.001 ***
#   Dilution        4 0.001620   3.6728  0.001 ***
#   OriginalSample 26 0.049452  17.2519  0.001 ***
#   Residual       71 0.007828                    


fig <-ordiplot(mod_apriori,type="n")
fig <-ordiplot(mod_apriori,type="n")
colvec <- c("black","grey51","blue","red","green") 


