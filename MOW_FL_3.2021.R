# load libraries ####
library(dismo)
library(gbm)
library(reshape2)
library(dplyr)
library(plyr)
library(xlsx)
library(ggplot2)
library(stringr)
library(ape)
library(scales)
library(data.table)

options(scipen=999)

# Install ggBRT 
install.packages("devtools") # if package "devtools" not already installed
devtools::install_github("JBjouffray/ggBRT")
library(ggBRT)

# set working directory and load data ####
setwd("")
rvc_sites <- read.csv("RVC_impact.csv", sep=",", quote="", na.strings=c("NA", "", " "), header=TRUE)

# separate impact model sites from biomass model sites ####

impact <- rvc_sites[rvc_sites$model=="impact",]
biomass <- rvc_sites[rvc_sites$model=="biomass",]


# collinearity ####
# We use pairwise relationship correlation coefficients and 
# variance inflation factor estimates to assess collinearity among predictors. 

# Select the set of numeric predictors
predictors <- c(27:71)
predictors.numeric <- c(27,29,30,31,33,34,35,36,37,40,41,42,43,44,45,46,48,
                        49,50,51,52,53,54,55,56,62,63,64,65,66,67,68,69,70) # all of the numeric predictors to test for colinearity

# Assess collinearity with Pearson correlation coefficient
library(corrplot)
coeff <- cor(impact[,predictors.numeric], method="pearson",use="pairwise.complete.obs")
coeff[which(abs(coeff[]) > 0.8)] # no correlation > 0.8

# Correlation plot
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(coeff,method="number", col=col(200),tl.cex=0.8,number.cex = 0.7,
         tl.col="black", tl.srt=45)

predictors.culled <- c(27,31,33,34,35,36,37,40,43,44,45,50,52,53,54,56,62,63,64,65,68)
coeff1 <- cor(impact[,predictors.culled], method="pearson",use="pairwise.complete.obs")
coeff1[which(abs(coeff1[]) > 0.8)] 
corrplot(coeff1,method="number", col=col(200),tl.cex=0.8,number.cex = 0.7,
         tl.col="black", tl.srt=45)


# Assess collinearity with Variable Inflation Factor (VIF)
#source("vif_func.R") # available at https://gist.github.com/fawda123/4717702#file-vif_fun-r
library(fmsb)

vif_func<-function(in_frame,thresh=10,trace=T,...){
  
  library(fmsb)
  
  if(any(!'data.frame' %in% class(in_frame))) in_frame<-data.frame(in_frame)
  
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  var_names <- names(in_frame)
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = in_frame, ...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]), na.rm = TRUE)
  
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(var_names)
  }
  else{
    
    in_dat<-in_frame
    
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      
      vif_vals<-NULL
      var_names <- names(in_dat)
      
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- formula(paste(val, '~', form))
        vif_add<-VIF(lm(form_in, data = in_dat, ...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2]), na.rm = TRUE))[1]
      
      vif_max<-as.numeric(vif_vals[max_row,2])
      
      if(vif_max<thresh) break
      
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
      
    }
    
    return(names(in_dat))
    
  }
  
}


vif_func(in_frame=rvc_sites[,predictors.culled],thresh=12,trace=T) # all variables have VIF < 12

rm(vif_func, col, predictors, predictors.culled, 
   predictors.numeric, coeff, coeff1)

# test regression tree parameters. RESULTS = use lr=0.01, bf=0.75, tc=5 ####

tree_complexity <- 1:5
learning_rate <- c(0.01,0.05,0.001,0.0001)
bag_fraction <- c(0.5,0.75,0.9)

names <- c("tc", "lr", "bf", "tot_dev", "resid_dev", "corr",
           "AUC", "perc_expl", "cv_dev", "cv_corr", "cv_AUC", "cv_perc_expl")

# data frame to store results
parameter_tests <- data.frame(matrix(NA, nrow = 1, ncol = 12))
colnames(parameter_tests) <- names

output <- data.frame(matrix(NA, nrow = 1, ncol = 12))
colnames(output) <- names

for (i in 1:length(tree_complexity)) {
  for (j in 1:length(learning_rate)) {
    for (k in 1:length(bag_fraction)) {
      
      model <- gbm.step(data = impact, 
                        gbm.x = c(27,31,33,34,35,36,37,40,43,44,45,50,52,53,54,56,62,63,64,65,68,71), 
                        gbm.y = 5,
                        family = "gaussian", 
                        tree.complexity = tree_complexity[i],
                        learning.rate = learning_rate[j],
                        bag.fraction = bag_fraction[k])
      
      output$tc <- tree_complexity[i]
      output$lr <- learning_rate[j]
      output$bf <- bag_fraction[k]
      output$tot_dev <- model$self.statistics$mean.null
      output$resid_dev <- model$self.statistics$mean.resid
      output$corr <- model$self.statistics$correlation
      output$AUC <- model$self.statistics$discrimination
      output$perc_expl <- (1 - model$self.statistics$mean.resid / model$self.statistics$mean.null)*100
      output$cv_dev <- model$cv.statistics$deviance.mean
      output$cv_corr <- model$cv.statistics$correlation.mean
      output$cv_AUC <- model$cv.statistics$discrimination.mean
      output$cv_perc_expl <- (1 - model$cv.statistics$deviance.mean / model$self.statistics$mean.null)*100
      
      parameter_tests <- rbind(parameter_tests, output)
      
    }
  }
}

rm(bag_fraction, i, j, k, learning_rate, names, tree_complexity)

# BRT impact model (response: SG biomass) ####

library(dismo)
library(ggBRT)
#https://github.com/JBjouffray/Hawaii_RegimesPredictors/blob/master/Scripts/Hawaii_RegimesPredictors.Rmd

impact.coral <- impact
impact.coral <- impact.coral[!impact.coral$Habitat_type_classLV2=="Pavement",]
impact.coral <- impact.coral[!impact.coral$Habitat_type_classLV2=="Pavement with Sand Channels",]
impact.coral <- impact.coral[!impact.coral$Habitat_type_classLV2=="Colonized Pavement",]
impact.coral <- impact.coral[!impact.coral$Habitat_type_classLV2=="Reef Rubble",]

### IMPACT MODEL (only coral reef habitats; excluding non-accreting hardbottom) ###
set.seed(4)
brt_impact.SGbiomass.coral <- gbm.step(data = impact.coral, 
                                       gbm.x = c(27,31,33,34,35,36,37,39,40,43,44,46,47,50,52,53,54,56,62,63,64,65,68,71), 
                                       gbm.y = 5, # log plus one (snapper-grouper biomass)
                                       family = "gaussian", 
                                       tree.complexity = 5,
                                       learning.rate = 0.01, 
                                       bag.fraction = 0.75) 

# extract the performance of the model
Perf <- ggPerformance(SG.biomass=brt_impact.SGbiomass.coral) 
round(Perf,2)

# Relative influence of predictors
ggInfluence(brt_impact.SGbiomass.coral,col.bar = "deepskyblue4", 
            show.signif=FALSE,
            main="Fishing impact model (Coral Reef only)")

# Partial dependency plots
ggPD(brt_impact.SGbiomass.coral, n.plots=12, col.line="black", cex.line=0.5, rug=TRUE,
     rug.pos="t", y.label="SG Biomass")
ggPDfit(brt_impact.SGbiomass.coral, n.plots=12, col.dot="grey", cex.dot=1)


### Significant Impact Model (Coral Reef habitats only)
set.seed(4)
brt_impact.SGbiomass.coral.signif <- gbm.step(data = impact.coral, 
                                              gbm.x = c(31,34,35,36,37,40,43,44,46,50,53,63), 
                                              gbm.y = 5, 
                                              family = "gaussian", 
                                              tree.complexity = 4,
                                              learning.rate = 0.01, 
                                              bag.fraction = 0.75)
Perf <- ggPerformance(SG.biomass=brt_impact.SGbiomass.coral.signif)
round(Perf,2)

ggInfluence(brt_impact.SGbiomass.coral.signif,col.bar = "deepskyblue4", 
            show.signif=FALSE,
            main="Fishing impact model (Coral Reef only)")

#Partial dependency plots (using ggBRT) 
brt_impact.SGbiomass.coral.signif.prerun<- plot.gbm.4list(brt_impact.SGbiomass.coral.signif)

# Boostrap the BRTs 1000 times to build confidence intervals
brt_impact.SGbiomass.coral.signif.boot <- gbm.bootstrap.functions(brt_impact.SGbiomass.coral.signif, 
          list.predictors=brt_impact.SGbiomass.coral.signif.prerun, n.reps=1000)

# Draw partial dependency plots for the x most influential predictors
ggPD_boot(brt_impact.SGbiomass.coral.signif,n.plots = 9, nrow = 1, list.4.preds = brt_impact.SGbiomass.coral.signif.prerun, 
          col.line="deepskyblue4", booted.preds = brt_impact.SGbiomass.coral.signif.boot$function.preds, 
          cex.line=1, col.ci="deepskyblue4", alpha.dot=0.2, type.ci = "ribbon", cex.axis=2, cis=c(0.05, 0.95),
          alpha.ci=0.1, rug = TRUE, y.label = "")

# Calculate Moran's I (spatial autocorrelation) for impact model
# convert lat/long values into a distance matrix (Euclidean distance)
# so that we know which sites are closest each other
library(geosphere)

distances.impact.coral <- as.matrix(dist(cbind(impact.coral$Latitude, impact.coral$Longitude)))
distances.impact.coral[1:5, 1:5]
distances.impact.coral.inverse <- 1/distances.impact.coral
diag(distances.impact.coral.inverse) <- 0
distances.impact.coral.inverse[1:5, 1:5]

# impact model residuals
head(brt_impact.SGbiomass.coral.signif$residuals)

# Moran.I(x=a numeric vector, weight=a matrix of weights)
Moran.I(brt_impact.SGbiomass.coral.signif$residuals, distances.impact.coral.inverse, na.rm=TRUE, alt="two.sided")
rm(distances.impact.coral, distances.impact.coral.inverse)

# Extrapolate fishing impact to entire FL reef tract ####

# load data (all reef pixels on the FL reef tract with associated data for
# explanatory variables that are significant in the fishing impact model)
ReefPoints_impact <- read.csv("ReefPoints_impact.csv", sep=",", quote="", na.strings=c("NA", "", " "), header=TRUE)

# Process data; set biophysical variables to their mean
ReefPoints_impact$Depth <- ReefPoints_impact$Depth*-1
ReefPoints_impact$Depth <- mean(ReefPoints_impact$Depth, na.rm=TRUE)

ReefPoints_impact$NPP <- mean(ReefPoints_impact$NPP, na.rm=TRUE) 
ReefPoints_impact$Wave_exposure <- mean(ReefPoints_impact$Wave_exposure, na.rm=TRUE)

ReefPoints_impact$Recreational_fishermen_50km <- ifelse(is.na(ReefPoints_impact$Recreational_fishermen_50km), 0, ReefPoints_impact$Recreational_fishermen_50km)

ReefPoints_impact$Coral_area_UFRTM_20km <- mean(ReefPoints_impact$Coral_area_UFRTM_20km, na.rm=TRUE)
ReefPoints_impact$Coral_area_UFRTM_20km <- as.integer(ReefPoints_impact$Coral_area_UFRTM_20km)

ReefPoints_impact$Deepwater <- mean(ReefPoints_impact$Deepwater, na.rm=TRUE)

ReefPoints_impact$Marina_slips_25km[is.na(ReefPoints_impact$Marina_slips_25km)] <- 0
ReefPoints_impact$Marina_slips_25km <- as.integer(ReefPoints_impact$Marina_slips_25km)

ReefPoints_impact$Nursery_mangroves <- mean(ReefPoints_impact$Nursery_mangroves, na.rm=TRUE)
ReefPoints_impact$SST <- ReefPoints_impact$SST - 273.15
ReefPoints_impact$SST <- mean(ReefPoints_impact$SST, na.rm=TRUE)

ReefPoints_impact$FSA <- mean(ReefPoints_impact$FSA, na.rm=TRUE)

ReefPoints_impact$Reef_complexity <- mean(impact.coral$Reef_complexity, na.rm=TRUE)
ReefPoints_impact$Coral_cover <- mean(impact.coral$Coral_cover, na.rm=TRUE)


# run predictions #
preds.impact <- predict.gbm(brt_impact.SGbiomass.coral.signif, ReefPoints_impact,
                            n.trees=brt_impact.SGbiomass.coral.signif$gbm.call$best.trees, 
                            type="response")

# rescale predictive values from 0 to 1 #
preds.impact <- rescale(preds.impact, to = c(1,0), from = range(preds.impact, na.rm = TRUE, finite = TRUE))

# amend prediction values back to pointid data table #
ReefPoints_impact$prediction.impact <- preds.impact
fishing_impact <- ReefPoints_impact[,c(1,16)]






# BRT biomass models ####

# load data #
# RVC sites for biomass model with fishing impact estimates (added to each rvc site in ArcGIS)
RVC_biomass <- read.csv("RVC_biomass.csv", sep=",", quote="", na.strings=c("NA", "", " "), header=TRUE)
biomass <- join(biomass, RVC_biomass, by="site") # adding the variable 'impact' to the biomass model dataset
rm(RVC_biomass)

# for Zuercher et al. 2021, use only data for coral reef habitats (not hardbottom habitats)
biomass <- biomass[!biomass$Habitat_type_classLV2=="Pavement",]
biomass <- biomass[!biomass$Habitat_type_classLV2=="Pavement with Sand Channels",]
biomass <- biomass[!biomass$Habitat_type_classLV2=="Colonized Pavement",]
biomass <- biomass[!biomass$Habitat_type_classLV2=="Reef Rubble",]

# SNAPPER-GROUPER BIOMASS #
set.seed(4)
brt_biomass.SG <- gbm.step(data = biomass, 
                           gbm.x = c(27,28,31,33,34,35,36,37,39,40,43,44,47,52,62,63,64,71,72), 
                           gbm.y = 5, 
                           family = "gaussian", 
                           tree.complexity = 5,
                           learning.rate = 0.01, 
                           bag.fraction = 0.75) 

Perf <- ggPerformance(SG.biomass=brt_biomass.SG) 
round(Perf,2)
ggInfluence(brt_biomass.SG,col.bar = "deepskyblue4", 
            show.signif=FALSE,
            main="Snapper-Grouper Current Biomass")

# run BRT model with just the important variables for SG biomass #
brt_biomass.SG.signif <- gbm.step(data = biomass, 
                                  gbm.x = c(28,31,34,35,36,37,40,43,63,64,72), 
                                  gbm.y = 5, 
                                  family = "gaussian", 
                                  tree.complexity = 5,
                                  learning.rate = 0.01, 
                                  bag.fraction = 0.75)
Perf <- ggPerformance(SG.biomass=brt_biomass.SG.signif) 
round(Perf,2)
ggInfluence(brt_biomass.SG.signif,col.bar = "deepskyblue4", 
            show.signif=FALSE,
            main="Snapper-Grouper Current Biomass")

brt_biomass.SG.signif.prerun<- plot.gbm.4list(brt_biomass.SG.signif)
brt_biomass.SG.signif.boot <- gbm.bootstrap.functions(brt_biomass.SG.signif, 
          list.predictors=brt_biomass.SG.signif.prerun, n.reps=1000)

ggPD_boot(brt_biomass.SG.signif,n.plots = 9, nrow = 1, list.4.preds = brt_biomass.SG.signif.prerun, 
          col.line="deepskyblue4", booted.preds = brt_biomass.SG.signif.boot$function.preds, 
          cex.line=1, col.ci="deepskyblue4", alpha.dot=0.2, type.ci = "ribbon", cex.axis=2, cis=c(0.05, 0.95),
          alpha.ci=0.1, rug = TRUE, y.label = "")

# Interactions between predictors #
ggInteract_list(brt_biomass.SG.signif) 

ggInteract_3D(gbm.object = brt_biomass.SG.signif, x="impact", y="Reef_complexity",col.gradient = "Blues",
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,
              x.range=c(0,1), y.range=c(0,12), z.range=c(5,7),
              label.contour = T,col.contour = "#92271F",show.axis = T,legend = T)



# Moran's I for SG biomass model #
distances <- as.matrix(dist(cbind(biomass$Latitude, biomass$Longitude)))
distances[1:5, 1:5]
distances.inverse <- 1/distances
diag(distances.inverse) <- 0
distances.inverse[1:5, 1:5]

# model residuals
head(brt_biomass.SG.signif$residuals)

# Moran's I for SG biomass model
Moran.I(brt_biomass.SG.signif$residuals, distances.inverse, na.rm=TRUE)


# TOTAL BIOMASS #
set.seed(4)
brt_biomass.total <- gbm.step(data = biomass, 
                              gbm.x = c(27,31,33,34,35,36,37,39,40,43,44,47,52,62,63,64,71,72), 
                              gbm.y = 3, 
                              family = "gaussian", 
                              tree.complexity = 5,
                              learning.rate = 0.01, 
                              bag.fraction = 0.75) 

Perf <- ggPerformance(total.biomass=brt_biomass.total) 
round(Perf,2)
ggInfluence(brt_biomass.total,col.bar = "deepskyblue4", 
            show.signif=FALSE,
            main="Total Current Biomass")

# run BRT model with just the important variables for total biomass #
set.seed(4)
brt_biomass.total.signif <- gbm.step(data = biomass, 
                                     gbm.x = c(31,34,36,37,40,43,72), 
                                     gbm.y = 3, 
                                     family = "gaussian", 
                                     tree.complexity = 5,
                                     learning.rate = 0.01, 
                                     bag.fraction = 0.75)
Perf <- ggPerformance(total.biomass=brt_biomass.total.signif) 
round(Perf,2)

ggInfluence(brt_biomass.total.signif,col.bar = "deepskyblue4", 
            show.signif=FALSE,
            main="Total Current Biomass")

brt_biomass.total.signif.prerun<- plot.gbm.4list(brt_biomass.total.signif)
brt_biomass.total.signif.boot <- gbm.bootstrap.functions(brt_biomass.total.signif, 
          list.predictors=brt_biomass.total.signif.prerun, n.reps=1000)

ggPD_boot(brt_biomass.total.signif, n.plots = 7, nrow = 1, list.4.preds = brt_biomass.total.signif.prerun, 
          col.line="deepskyblue4", booted.preds = brt_biomass.total.signif.boot$function.preds, 
          cex.line=1, col.ci="deepskyblue4", alpha.dot=0.2, type.ci = "ribbon", cex.axis=2, cis=c(0.05, 0.95),
          alpha.ci=0.1, rug = TRUE, y.label = "")


# Moran's I for total biomass model
Moran.I(brt_biomass.total.signif$residuals, distances.inverse, na.rm=TRUE)


# MARINE LIFE BIOMASS
set.seed(10)
brt_biomass.ML <- gbm.step(data = biomass, 
                           gbm.x = c(27,31,33,34,35,36,37,39,40,43,44,47,52,62,63,64,71,72), 
                           gbm.y = 7, 
                           family = "gaussian", 
                           tree.complexity = 5,
                           learning.rate = 0.01, 
                           bag.fraction = 0.75) 

Perf <- ggPerformance(ML.biomass=brt_biomass.ML) 
round(Perf,2)
ggInfluence(brt_biomass.ML,col.bar = "deepskyblue4", 
            show.signif=FALSE,
            main="Marine Life Species Current Biomass")

# run BRT model with just the important variables for marine life species biomass #
set.seed(10)
brt_biomass.ML.signif <- gbm.step(data = biomass, 
                                  gbm.x = c(27,31,33,34,35,36,43,72), 
                                  gbm.y = 7, 
                                  family = "gaussian", 
                                  tree.complexity = 5,
                                  learning.rate = 0.01, 
                                  bag.fraction = 0.75)
Perf <- ggPerformance(ML.biomass=brt_biomass.ML.signif) 
round(Perf,2)

ggInfluence(brt_biomass.ML.signif,col.bar = "deepskyblue4", 
            show.signif=FALSE,
            main="Marine Life Current Biomass")

brt_biomass.ML.signif.prerun<- plot.gbm.4list(brt_biomass.ML.signif)
brt_biomass.ML.signif.boot <- gbm.bootstrap.functions(brt_biomass.ML.signif, 
          list.predictors=brt_biomass.ML.signif.prerun, n.reps=1000)

ggPD_boot(brt_biomass.ML.signif,n.plots = 8, nrow = 1, list.4.preds = brt_biomass.ML.signif.prerun, 
          col.line="deepskyblue4", booted.preds = brt_biomass.ML.signif.boot$function.preds, 
          cex.line=1, col.ci="deepskyblue4", alpha.dot=0.2, type.ci = "ribbon", cex.axis=2, cis=c(0.05, 0.95),
          alpha.ci=0.1, rug = TRUE, y.label = "")


# Moran's I for marine life biomass model
Moran.I(brt_biomass.ML.signif$residuals, distances.inverse, na.rm=TRUE)


# PARROTFISH BIOMASS
set.seed(4)
brt_biomass.pf <- gbm.step(data = biomass, 
                           gbm.x = c(27,31,33,34,35,36,37,39,40,43,44,47,52,62,63,64,71,72), 
                           gbm.y = 8, 
                           family = "gaussian", 
                           tree.complexity = 5,
                           learning.rate = 0.01, 
                           bag.fraction = 0.75) 

Perf <- ggPerformance(pf.biomass=brt_biomass.pf) 
round(Perf,2)
ggInfluence(brt_biomass.pf,col.bar = "deepskyblue4", 
            show.signif=FALSE,
            main="Parrotfish Current Biomass")

# run BRT model with just the important variables for parrotfish biomass #
set.seed(4)
brt_biomass.pf.signif <- gbm.step(data = biomass, 
                                  gbm.x = c(27,31,33,34,35,36,37,39,40,43,62,72), 
                                  gbm.y = 8, 
                                  family = "gaussian", 
                                  tree.complexity = 5,
                                  learning.rate = 0.01, 
                                  bag.fraction = 0.75)
Perf <- ggPerformance(pf.biomass=brt_biomass.pf.signif) 
round(Perf,2)

ggInfluence(brt_biomass.pf.signif,col.bar = "deepskyblue4", 
            show.signif=FALSE,
            main="Parrotfish Current Biomass")

brt_biomass.pf.signif.prerun<- plot.gbm.4list(brt_biomass.pf.signif)
brt_biomass.pf.signif.boot <- gbm.bootstrap.functions(brt_biomass.pf.signif, 
          list.predictors=brt_biomass.pf.signif.prerun, n.reps=1000)

ggPD_boot(brt_biomass.pf.signif,n.plots = 8, nrow = 1, list.4.preds = brt_biomass.pf.signif.prerun, 
          col.line="deepskyblue4", booted.preds = brt_biomass.pf.signif.boot$function.preds, 
          cex.line=1, col.ci="deepskyblue4", alpha.dot=0.2, type.ci = "ribbon", cex.axis=2, cis=c(0.05, 0.95),
          alpha.ci=0.1, rug = TRUE, y.label = "")

# Moran's I for parrotfish biomass model
Moran.I(brt_biomass.pf.signif$residuals, distances.inverse, na.rm=TRUE)



# HERBIVORE BIOMASS
set.seed(4)
brt_biomass.herbi <- gbm.step(data = biomass, 
                              gbm.x = c(27,31,33,34,35,36,37,39,40,43,44,47,52,62,63,64,71,72), 
                              gbm.y = 6, 
                              family = "gaussian", 
                              tree.complexity = 5,
                              learning.rate = 0.01, 
                              bag.fraction = 0.75) 

Perf <- ggPerformance(herbi.biomass=brt_biomass.herbi) 
round(Perf,2)
ggInfluence(brt_biomass.herbi,col.bar = "deepskyblue4", 
            show.signif=FALSE,
            main="Herbivore Current Biomass")

# run BRT model with just the important variables for parrotfish biomass #
set.seed(4)
brt_biomass.herbi.signif <- gbm.step(data = biomass, 
                                     gbm.x = c(27,31,33,34,35,36,37,39,43,72), 
                                     gbm.y = 6, 
                                     family = "gaussian", 
                                     tree.complexity = 5,
                                     learning.rate = 0.01, 
                                     bag.fraction = 0.75)
Perf <- ggPerformance(herbi.biomass=brt_biomass.herbi.signif) 
round(Perf,2)

ggInfluence(brt_biomass.herbi.signif,col.bar = "deepskyblue4", 
            show.signif=FALSE,
            main="Herbivore Current Biomass")

brt_biomass.herbi.signif.prerun<- plot.gbm.4list(brt_biomass.herbi.signif)
brt_biomass.herbi.signif.boot <- gbm.bootstrap.functions(brt_biomass.herbi.signif, 
          list.predictors=brt_biomass.herbi.signif.prerun, n.reps=1000)

ggPD_boot(brt_biomass.herbi.signif,n.plots = 10, nrow = 1, list.4.preds = brt_biomass.herbi.signif.prerun, 
          col.line="deepskyblue4", booted.preds = brt_biomass.herbi.signif.boot$function.preds, 
          cex.line=1, col.ci="deepskyblue4", alpha.dot=0.2, type.ci = "ribbon", cex.axis=2, cis=c(0.05, 0.95),
          alpha.ci=0.1, rug = TRUE, y.label = "")

# Moran's I for parrotfish biomass model
Moran.I(brt_biomass.herbi.signif$residuals, distances.inverse, na.rm=TRUE)





# Extrapolate current biomass to entire FL reef tract ####

# data frame with all of the Coral Reef grid points (i.e. points at the center of 
# each reef pixel/cell) with appended attributes for each important variable in the biomass model
ReefPoints_biomass <- read.csv("ReefPoints_biomass.csv", sep=",", quote="", na.strings=c("NA", "", " "), header=TRUE)

# Dataframe for Current Biomass results
Current_biomass <- ReefPoints_biomass

# Predictive BRT: SG BIOMASS #
preds.SG <- predict.gbm(brt_biomass.SG.signif, Current_biomass,
                        n.trees=brt_biomass.SG.signif$gbm.call$best.trees, 
                        type="response")

# amend biomass prediction values back to pointid data table #
Current_biomass$current_biomass.SG <- preds.SG
# revert back to actual biomass values rather than logp1
Current_biomass$current_biomass.SG <- exp(Current_biomass$current_biomass.SG) - 1


# Predictive BRT: TOTAL BIOMASS #
preds.total <- predict.gbm(brt_biomass.total.signif, Current_biomass,
                           n.trees=brt_biomass.total.signif$gbm.call$best.trees, 
                           type="response")
Current_biomass$current_biomass.total <- preds.total
Current_biomass$current_biomass.total <- exp(Current_biomass$current_biomass.total) - 1


# Predictive BRT: PARROTFISH BIOMASS #
preds.pf <- predict.gbm(brt_biomass.pf.signif, Current_biomass,
                        n.trees=brt_biomass.pf.signif$gbm.call$best.trees, 
                        type="response")
Current_biomass$current_biomass.pf <- preds.pf
Current_biomass$current_biomass.pf <- exp(Current_biomass$current_biomass.pf) - 1


# Predictive BRT: HERBI BIOMASS #
preds.herbi <- predict.gbm(brt_biomass.herbi.signif, Current_biomass,
                           n.trees=brt_biomass.herbi.signif$gbm.call$best.trees, 
                           type="response")
Current_biomass$current_biomass.herbi <- preds.herbi
Current_biomass$current_biomass.herbi <- exp(Current_biomass$current_biomass.herbi) - 1



# Predictive BRT: MARINE LIFE SPECIES BIOMASS #
preds.ML <- predict.gbm(brt_biomass.ML.signif, Current_biomass,
                        n.trees=brt_biomass.ML.signif$gbm.call$best.trees, 
                        type="response")
Current_biomass$current_biomass.ML <- preds.ML
Current_biomass$current_biomass.ML <- exp(Current_biomass$current_biomass.ML) - 1





# Reef restoration scenarios ####

# SCENARIO 1
# NOAA Reef Restoration (FL Iconic Reefs) Phase 1a
Scenario1 <- ReefPoints_biomass
Scenario1$Coral_cover <- 10
Scenario1$Reef_complexity <- Scenario1$Reef_complexity + 0.15

preds.S1.SG <- predict.gbm(brt_biomass.SG.signif, Scenario1,
                           n.trees=brt_biomass.SG.signif$gbm.call$best.trees, type="response")
Scenario1$S1.SG <- preds.S1.SG
Scenario1$S1.SG <- exp(Scenario1$S1.SG) - 1

preds.S1.total <- predict.gbm(brt_biomass.total.signif, Scenario1,
                              n.trees=brt_biomass.SG.signif$gbm.call$best.trees, type="response")
Scenario1$S1.total <- preds.S1.total
Scenario1$S1.total <- exp(Scenario1$S1.total) - 1

preds.S1.herbi <- predict.gbm(brt_biomass.herbi.signif, Scenario1,
                              n.trees=brt_biomass.herbi.signif$gbm.call$best.trees, type="response")
Scenario1$S1.herbi <- preds.S1.herbi
Scenario1$S1.herbi <- exp(Scenario1$S1.herbi) - 1

preds.S1.ML <- predict.gbm(brt_biomass.ML.signif, Scenario1,
                           n.trees=brt_biomass.ML.signif$gbm.call$best.trees, type="response")
Scenario1$S1.ML <- preds.S1.ML
Scenario1$S1.ML <- exp(Scenario1$S1.ML) - 1



# SCENARIO 2
# NOAA Reef Restoration (FL Iconic Reefs) Phase 2
Scenario2 <- ReefPoints_biomass
Scenario2$Coral_cover <- 25
Scenario2$Reef_complexity <- Scenario2$Reef_complexity + 0.75

preds.S2.SG <- predict.gbm(brt_biomass.SG.signif, Scenario2,
                           n.trees=brt_biomass.SG.signif$gbm.call$best.trees, type="response")
Scenario2$S2.SG <- preds.S2.SG
Scenario2$S2.SG <- exp(Scenario2$S2.SG) - 1

preds.S2.total <- predict.gbm(brt_biomass.total.signif, Scenario2,
                              n.trees=brt_biomass.SG.signif$gbm.call$best.trees, type="response")
Scenario2$S2.total <- preds.S2.total
Scenario2$S2.total <- exp(Scenario2$S2.total) - 1

preds.S2.herbi <- predict.gbm(brt_biomass.herbi.signif, Scenario2,
                              n.trees=brt_biomass.herbi.signif$gbm.call$best.trees, type="response")
Scenario2$S2.herbi <- preds.S2.herbi
Scenario2$S2.herbi <- exp(Scenario2$S2.herbi) - 1

preds.S2.ML <- predict.gbm(brt_biomass.ML.signif, Scenario2,
                           n.trees=brt_biomass.ML.signif$gbm.call$best.trees, type="response")
Scenario2$S2.ML <- preds.S2.ML
Scenario2$S2.ML <- exp(Scenario2$S2.ML) - 1




# SCENARIO 3
# Addition of artificial reef structure
Scenario3 <- ReefPoints_biomass
Scenario3$Reef_complexity <- Scenario3$Reef_complexity + 1.2

preds.S3.SG <- predict.gbm(brt_biomass.SG.signif, Scenario3,
                           n.trees=brt_biomass.SG.signif$gbm.call$best.trees, type="response")
Scenario3$S3.SG <- preds.S3.SG
Scenario3$S3.SG <- exp(Scenario3$S3.SG) - 1

preds.S3.total <- predict.gbm(brt_biomass.total.signif, Scenario3,
                              n.trees=brt_biomass.SG.signif$gbm.call$best.trees, type="response")
Scenario3$S3.total <- preds.S3.total
Scenario3$S3.total <- exp(Scenario3$S3.total) - 1

preds.S3.herbi <- predict.gbm(brt_biomass.herbi.signif, Scenario3,
                              n.trees=brt_biomass.herbi.signif$gbm.call$best.trees, type="response")
Scenario3$S3.herbi <- preds.S3.herbi
Scenario3$S3.herbi <- exp(Scenario3$S3.herbi) - 1

preds.S3.ML <- predict.gbm(brt_biomass.ML.signif, Scenario3,
                           n.trees=brt_biomass.ML.signif$gbm.call$best.trees, type="response")
Scenario3$S3.ML <- preds.S3.ML
Scenario3$S3.ML <- exp(Scenario3$S3.ML) - 1





# SCENARIO 4
# No fishing
Scenario4 <- ReefPoints_biomass
Scenario4$impact <- 0

preds.S4.SG <- predict.gbm(brt_biomass.SG.signif, Scenario4,
                           n.trees=brt_biomass.SG.signif$gbm.call$best.trees, type="response")
Scenario4$S4.SG <- preds.S4.SG
Scenario4$S4.SG <- exp(Scenario4$S4.SG) - 1

preds.S4.total <- predict.gbm(brt_biomass.total.signif, Scenario4,
                              n.trees=brt_biomass.SG.signif$gbm.call$best.trees, type="response")
Scenario4$S4.total <- preds.S4.total
Scenario4$S4.total <- exp(Scenario4$S4.total) - 1

preds.S4.herbi <- predict.gbm(brt_biomass.herbi.signif, Scenario4,
                              n.trees=brt_biomass.herbi.signif$gbm.call$best.trees, type="response")
Scenario4$S4.herbi <- preds.S4.herbi
Scenario4$S4.herbi <- exp(Scenario4$S4.herbi) - 1

preds.S4.ML <- predict.gbm(brt_biomass.ML.signif, Scenario4,
                           n.trees=brt_biomass.ML.signif$gbm.call$best.trees, type="response")
Scenario4$S4.ML <- preds.S4.ML
Scenario4$S4.ML <- exp(Scenario4$S4.ML) - 1






# SCENARIO 5
# No Fishing + Reef Restoration Phase 2
Scenario5 <- ReefPoints_biomass
Scenario5$Coral_cover <- 25
Scenario5$impact <- 0
Scenario5$Reef_complexity <- Scenario5$Reef_complexity + 0.75

preds.S5.SG <- predict.gbm(brt_biomass.SG.signif, Scenario5,
                           n.trees=brt_biomass.SG.signif$gbm.call$best.trees, type="response")
Scenario5$S5.SG <- preds.S5.SG
Scenario5$S5.SG <- exp(Scenario5$S5.SG) - 1

preds.S5.total <- predict.gbm(brt_biomass.total.signif, Scenario5,
                              n.trees=brt_biomass.SG.signif$gbm.call$best.trees, type="response")
Scenario5$S5.total <- preds.S5.total
Scenario5$S5.total <- exp(Scenario5$S5.total) - 1

preds.S5.herbi <- predict.gbm(brt_biomass.herbi.signif, Scenario5,
                              n.trees=brt_biomass.herbi.signif$gbm.call$best.trees, type="response")
Scenario5$S5.herbi <- preds.S5.herbi
Scenario5$S5.herbi <- exp(Scenario5$S5.herbi) - 1

preds.S5.ML <- predict.gbm(brt_biomass.ML.signif, Scenario5,
                           n.trees=brt_biomass.ML.signif$gbm.call$best.trees, type="response")
Scenario5$S5.ML <- preds.S5.ML
Scenario5$S5.ML <- exp(Scenario5$S5.ML) - 1




# SCENARIO 6
# No Fishing + Artificial structure (III)
Scenario6 <- ReefPoints_biomass
Scenario6$impact <- 0
Scenario6$Reef_complexity <- Scenario6$Reef_complexity + 1.2

preds.S6.SG <- predict.gbm(brt_biomass.SG.signif, Scenario6,
                           n.trees=brt_biomass.SG.signif$gbm.call$best.trees, type="response")
Scenario6$S6.SG <- preds.S6.SG
Scenario6$S6.SG <- exp(Scenario6$S6.SG) - 1

preds.S6.total <- predict.gbm(brt_biomass.total.signif, Scenario6,
                              n.trees=brt_biomass.SG.signif$gbm.call$best.trees, type="response")
Scenario6$S6.total <- preds.S6.total
Scenario6$S6.total <- exp(Scenario6$S6.total) - 1

preds.S6.herbi <- predict.gbm(brt_biomass.herbi.signif, Scenario6,
                              n.trees=brt_biomass.herbi.signif$gbm.call$best.trees, type="response")
Scenario6$S6.herbi <- preds.S6.herbi
Scenario6$S6.herbi <- exp(Scenario6$S6.herbi) - 1

preds.S6.ML <- predict.gbm(brt_biomass.ML.signif, Scenario6,
                           n.trees=brt_biomass.ML.signif$gbm.call$best.trees, type="response")
Scenario6$S6.ML <- preds.S6.ML
Scenario6$S6.ML <- exp(Scenario6$S6.ML) - 1




# Final results dataset ####

Results <- Current_biomass[,c(1,10:14,24:28)]
Results <- join(Results, Scenario1[,c(1,24:27)], by="pointid")
Results <- join(Results, Scenario2[,c(1,24:27)], by="pointid")
Results <- join(Results, Scenario3[,c(1,24:27)], by="pointid")
Results <- join(Results, Scenario4[,c(1,24:27)], by="pointid")
Results <- join(Results, Scenario5[,c(1,24:27)], by="pointid")
Results <- join(Results, Scenario6[,c(1,24:27)], by="pointid")

# to calculate results by region
BNP <- Results[!is.na(Results$BNP),]
CoralECA <- Results[!is.na(Results$CoralECA),]
FKNMS <- Results[!is.na(Results$FKNMS),]
DryTortugas <- Results[!is.na(Results$DryTortugas),]

# calculate mean and sd for each regional scenario

mean <- apply(CoralECA[,c(7:35)], 2, mean)
sd <- apply(CoralECA[,c(7:35)], 2, sd)
CoralECA_summary <- data.frame(mean, sd)

mean <- apply(BNP[,c(7:35)], 2, mean)
sd <- apply(BNP[,c(7:35)], 2, sd)
BNP_summary <- data.frame(mean, sd)

mean <- apply(FKNMS[,c(7:35)], 2, mean)
sd <- apply(FKNMS[,c(7:35)], 2, sd)
FKNMS_summary <- data.frame(mean, sd)

mean <- apply(DryTortugas[,c(7:35)], 2, mean)
sd <- apply(DryTortugas[,c(7:35)], 2, sd)
DryTortugas_summary <- data.frame(mean, sd)


# Export data layers for Scenario maps

# export #
fwrite(Results[,c(1,7,12,16,20,24,28,32)], "scenarios_SG.csv")



