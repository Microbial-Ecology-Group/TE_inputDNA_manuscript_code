### Raw reads ####
AMR_mapfile %>%
  summarise(total_raw_reads = sum(Raw_reads, na.rm = TRUE),
            mean_raw_reads = mean(Raw_reads, na.rm = TRUE),
            min_raw_reads = min(Raw_reads, na.rm = TRUE),
            max_raw_reads = max(Raw_reads, na.rm = TRUE))

AMR_mapfile <- AMR_mapfile %>%
  dplyr::mutate(AMR_percentage = (AMR_hits / Raw_reads) * 100)


AMR_mapfile %>%
  summarise(total_AMR_percentage = sum(AMR_percentage, na.rm = TRUE),
            mean_AMR_percentage = mean(AMR_percentage, na.rm = TRUE),
            min_AMR_percentage = min(AMR_percentage, na.rm = TRUE),
            max_AMR_percentage = max(AMR_percentage, na.rm = TRUE))



# Kit
pairwise.wilcox.test(AMR_mapfile$Raw_reads ,AMR_mapfile$Kit, p.adjust.method = "BH")

# InputDNA
pairwise.wilcox.test(AMR_mapfile$Raw_reads ,AMR_mapfile$InputDNA, p.adjust.method = "BH")

# SampleGroup
pairwise.wilcox.test(AMR_mapfile$Raw_reads ,AMR_mapfile$SampleGroup, p.adjust.method = "BH")
boxplot(AMR_mapfile$Raw_reads  ~AMR_mapfile$SampleGroup)


# OriginalSample
pairwise.wilcox.test(AMR_mapfile$Raw_reads ,AMR_mapfile$OriginalSample, p.adjust.method = "BH")





### AMR Hits ####
#### AMR hit sums #####
# total hits
sum(sample_sums(data_AMR.ps)) # 97334772

# average
mean(sample_sums(data_AMR.ps)) # 884861.6

# min, max
min(sample_sums(data_AMR.ps)) #127467
max(sample_sums(data_AMR.ps)) # 1379975

# Percent of hits to SNP confirmed reads
sum(sample_sums(data_SNPconfirm)) / sum(sample_sums(data_noSNP)) * 100


DatasetA_metadata <- subset(AMR_mapfile, SampleGroup != "Canine fecal")

# Kit
pairwise.wilcox.test(DatasetA_metadata$AMR_hits,DatasetA_metadata$Kit, p.adjust.method = "BH")
boxplot(DatasetA_metadata$AMR_hits ~DatasetA_metadata$Kit)

# InputDNA
pairwise.wilcox.test(DatasetA_metadata$AMR_hits,DatasetA_metadata$InputDNA, p.adjust.method = "BH")

# SampleGroup
pairwise.wilcox.test(DatasetA_metadata$AMR_hits,DatasetA_metadata$SampleGroup, p.adjust.method = "BH")
boxplot(DatasetA_metadata$AMR_hits ~DatasetA_metadata$SampleGroup)


# OriginalSample
pairwise.wilcox.test(DatasetA_metadata$AMR_hits,DatasetA_metadata$OriginalSample, p.adjust.method = "BH")


### By Species group, comparing Kit ####

#### Raw_reads ####
PorcineFecal_AMR_mapfile <- subset(AMR_mapfile, SampleGroup == "Porcine fecal")
pairwise.wilcox.test(PorcineFecal_AMR_mapfile$Raw_reads,PorcineFecal_AMR_mapfile$Kit, p.adjust.method = "BH")
boxplot(PorcineFecal_AMR_mapfile$Raw_reads ~PorcineFecal_AMR_mapfile$Kit)
# p 0.039

HumanFecal_AMR_mapfile <- subset(AMR_mapfile, SampleGroup == "Human wastewater")
pairwise.wilcox.test(HumanFecal_AMR_mapfile$Raw_reads,HumanFecal_AMR_mapfile$Kit, p.adjust.method = "BH")
boxplot(HumanFecal_AMR_mapfile$Raw_reads ~HumanFecal_AMR_mapfile$Kit)
# 5e-04

BovineFecal_AMR_mapfile <- subset(AMR_mapfile, SampleGroup == "Bovine fecal")
pairwise.wilcox.test(BovineFecal_AMR_mapfile$Raw_reads,BovineFecal_AMR_mapfile$Kit, p.adjust.method = "BH")
boxplot(HumanFecal_AMR_mapfile$Raw_reads ~HumanFecal_AMR_mapfile$Kit)


AvianFecal_AMR_mapfile <- subset(AMR_mapfile, SampleGroup == "Avian fecal")
pairwise.wilcox.test(AvianFecal_AMR_mapfile$Raw_reads,AvianFecal_AMR_mapfile$Kit, p.adjust.method = "BH")
boxplot(AvianFecal_AMR_mapfile$Raw_reads ~AvianFecal_AMR_mapfile$Kit)
# 0.0086



#### AMR_hits ####
PorcineFecal_AMR_mapfile <- subset(AMR_mapfile, SampleGroup == "Porcine fecal")
pairwise.wilcox.test(PorcineFecal_AMR_mapfile$AMR_hits,PorcineFecal_AMR_mapfile$Kit, p.adjust.method = "BH")

HumanFecal_AMR_mapfile <- subset(AMR_mapfile, SampleGroup == "Human wastewater")
pairwise.wilcox.test(HumanFecal_AMR_mapfile$AMR_hits,HumanFecal_AMR_mapfile$Kit, p.adjust.method = "BH")
#boxplot(HumanFecal_AMR_mapfile$AMR_hits ~HumanFecal_AMR_mapfile$Kit)

BovineFecal_AMR_mapfile <- subset(AMR_mapfile, SampleGroup == "Bovine fecal")
pairwise.wilcox.test(BovineFecal_AMR_mapfile$AMR_hits,BovineFecal_AMR_mapfile$Kit, p.adjust.method = "BH")

AvianFecal_AMR_mapfile <- subset(AMR_mapfile, SampleGroup == "Avian fecal")
pairwise.wilcox.test(AvianFecal_AMR_mapfile$AMR_hits,AvianFecal_AMR_mapfile$Kit, p.adjust.method = "BH")
# 0.023
boxplot(AvianFecal_AMR_mapfile$AMR_hits ~AvianFecal_AMR_mapfile$Kit)



### By Species group, comparing InputDNA ####

#### AMR_hits ####
CanineFecal_AMR_mapfile <- subset(AMR_mapfile, SampleGroup == "Canine fecal")
pairwise.wilcox.test(CanineFecal_AMR_mapfile$AMR_hits,CanineFecal_AMR_mapfile$InputDNA, p.adjust.method = "BH")

PorcineFecal_AMR_mapfile <- subset(AMR_mapfile, SampleGroup == "Porcine fecal")
pairwise.wilcox.test(PorcineFecal_AMR_mapfile$AMR_hits,PorcineFecal_AMR_mapfile$InputDNA, p.adjust.method = "BH")

HumanFecal_AMR_mapfile <- subset(AMR_mapfile, SampleGroup == "Human wastewater")
pairwise.wilcox.test(HumanFecal_AMR_mapfile$AMR_hits,HumanFecal_AMR_mapfile$InputDNA, p.adjust.method = "BH")
boxplot(HumanFecal_AMR_mapfile$AMR_hits ~HumanFecal_AMR_mapfile$InputDNA)
# 200 t0 800ng

BovineFecal_AMR_mapfile <- subset(AMR_mapfile, SampleGroup == "Bovine fecal")
pairwise.wilcox.test(BovineFecal_AMR_mapfile$AMR_hits,BovineFecal_AMR_mapfile$InputDNA, p.adjust.method = "BH")

AvianFecal_AMR_mapfile <- subset(AMR_mapfile, SampleGroup == "Avian fecal")
pairwise.wilcox.test(AvianFecal_AMR_mapfile$AMR_hits,AvianFecal_AMR_mapfile$InputDNA, p.adjust.method = "BH")

#### Raw_reads ####
CanineFecal_AMR_mapfile <- subset(AMR_mapfile, SampleGroup == "Canine fecal")
pairwise.wilcox.test(CanineFecal_AMR_mapfile$Raw_reads,CanineFecal_AMR_mapfile$InputDNA, p.adjust.method = "BH")

PorcineFecal_AMR_mapfile <- subset(AMR_mapfile, SampleGroup == "Porcine fecal")
pairwise.wilcox.test(PorcineFecal_AMR_mapfile$Raw_reads,PorcineFecal_AMR_mapfile$InputDNA, p.adjust.method = "BH")

HumanFecal_AMR_mapfile <- subset(AMR_mapfile, SampleGroup == "Human wastewater")
pairwise.wilcox.test(HumanFecal_AMR_mapfile$Raw_reads,HumanFecal_AMR_mapfile$InputDNA, p.adjust.method = "BH")

BovineFecal_AMR_mapfile <- subset(AMR_mapfile, SampleGroup == "Bovine fecal")
pairwise.wilcox.test(BovineFecal_AMR_mapfile$Raw_reads,BovineFecal_AMR_mapfile$InputDNA, p.adjust.method = "BH")

AvianFecal_AMR_mapfile <- subset(AMR_mapfile, SampleGroup == "Avian fecal")
pairwise.wilcox.test(AvianFecal_AMR_mapfile$Raw_reads,AvianFecal_AMR_mapfile$InputDNA, p.adjust.method = "BH")





# Run a generalized linear model with all variables #####
model <- glm(AMR_hits ~ Kit + Dilution + InputDNA + SampleGroup + Scaled_Raw_reads, data = AMR_mapfile)

# Print the summary of the model
summary(model)
model <- glm(AMR_hits ~ Kit + Dilution + InputDNA + SampleGroup, data = AMR_mapfile, family = poisson())

library(MASS)
library(AER)
nb_model <- glm.nb(AMR_hits ~ Kit + Dilution + InputDNA + SampleGroup + OriginalSample, data = AMR_mapfile)
confint(model)
summary(model)
# Fit the Poisson GLM
poisson_model <- glm(AMR_hits ~ Kit + Dilution + InputDNA + SampleGroup + OriginalSample, 
                     data = AMR_mapfile, family = poisson())

# Check for overdispersion
dispersion_test <- dispersiontest(poisson_model)

# Print the results
dispersion_test

#### generalized linear mixed-effects models ####
library(lme4)


# Fit a negative binomial mixed-effects model
library(glmmTMB)
library(MuMIn)
# Print the summary of the model

AMR_mapfile <- na.omit(AMR_mapfile)

# Fit the full model with glmer()
model_full <- glmmTMB(AMR_hits ~ Kit + Dilution + InputDNA + Scaled_Raw_reads + SampleGroup + (1 | OriginalSample), 
                      data = AMR_mapfile, family = nbinom2(), na.action = na.exclude)


model_reduced <- glmmTMB(AMR_hits ~ Kit  + InputDNA + Scaled_Raw_reads + SampleGroup + (1 | OriginalSample), 
                         data = AMR_mapfile, family = nbinom2(), na.action = na.exclude)


anova(model_reduced, model_full)

summary(model_full)


# Specify the model
model_full <- glmer(AMR_hits ~ Kit + Dilution + InputDNA  + (1|OriginalSample), 
                    data = AMR_mapfile, family = poisson())

# Reduced model without Kit
model_reduced <- lmer(AMR_hits ~ Dilution + InputDNA + (1|OriginalSample), 
                      data = AMR_mapfile, family = poisson())

anova(model_reduced, model_full)



# Check overdispersion
overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) sum(1:(dim(m)[1]))
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nobs(model)-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(model)

# Check random effects distribution
ranef_hist <- hist(ranef(model)$SampleGroup)

# Check functional form with residual plots
plot(resid(model)~fitted(model))






###### ANOVA #####
# Load necessary library
library(stats)
library(MASS)

# Remove rows with "1:10 dilution" in the "Dilution" column
filtered_data <- subset(AMR_mapfile, Dilution != "1:10 dilution")

# AMR_hits - Run ANOVA with multiple independent variables ####
model <- aov(AMR_hits ~ Kit + Dilution + InputDNA + SampleGroup + OriginalSample, data = filtered_data)

# Fit your initial ANOVA model
initial_model <- aov(AMR_hits ~ Kit + Dilution + InputDNA + SampleGroup + OriginalSample, data = filtered_data)

# Perform stepwise selection based on AIC
final_model <- stepAIC(initial_model, direction = "both")

# Print the final model summary
summary(final_model)

# Print the summary of the model
print(summary(final_model))

### Testing normality of residuals ####
# Extract the residuals
residuals <- residuals(final_model)

# Perform the Shapiro-Wilk test
shapiro_test <- shapiro.test(residuals)

# Print the test results, non sig so data does not deviate significantly from a normal distribution.
print(shapiro_test)

#### Homogeneity  of variance ####
# Extract the residuals from the ANOVA model
residuals <- residuals(final_model)

# Run the Levene's test on the residuals
# Dilution was significant (0.04816)
leveneTest(residuals, group = filtered_data$Dilution)


### eta squared
# Load the required package
library(effectsize)

# Calculate eta-squared
eta_sq <- eta_squared(final_model)

# Print the results
print(eta_sq)

#### Tukeys ####
# Perform pairwise comparisons using Tukey's HSD test
posthoc <- TukeyHSD(final_model)

# Print the results
print(posthoc)


### Plotting raw reads
ggplot(alpha_div_meta[alpha_div_meta$SampleGroup != "Canine fecal",], aes(x= SampleGroup, y= AMR_hits, fill = Kit, colour = Kit)) +
  theme_bw() + 
  labs(y= "AMR hits") +
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

