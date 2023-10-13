### load the libraries
# modified from Nuri's file of the same name
# ALP 3/12/2023
rm()

library(lme4)
library(emmeans)
library(sjPlot)
library(report)
library(ggplot2)


  #load MATLAB output (in table format in .txt) file in R
  foldN = "//ad.gatech.edu/bme/labs/singer/Abby/code/chronicflicker-ephys-prospectivecoding/results/LMM_R/"
  
 
  ##### BEHAVIOR - TRIAL #####
  ## load data file
  fileN = "TableData_Behavior_Trial.txt"
  dt = "Behavior_Trial"
  behavTypes <-c("duration", "lickDI", "nLicks", "licklatency_s", "licklatency_pos", "trial_speed", "AZ_speed", "RZ_speed", "Ctrl_speed", "AZ_speed_slope", "AZ_lick_slope")
  
  #initialize table structure to add stats to
  outputStats = data.frame()
  
  tblFilename = paste0(foldN, fileN)
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
  csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table
  
  # initialize text file
  close(file(txtFilename, open="w" ) )
  
  #read data
  mydata <- read.csv(tblFilename, head=T)

  #factorize variables to be treated as categories/levels instead of numerical value
  mydata$day = as.factor(mydata$day)
  mydata$animal = as.factor(mydata$animal)
  mydata$group = as.factor(mydata$group)

  ## get data subset of interest
  # only use rewarded, engaged, full trials
  testData <- mydata[mydata$engaged == 1 & #engaged
                         mydata$fullTrial == 1 & #full pass through the trial
                         mydata$rewarded == 1 & #received rewards on this trial
                         mydata$timepoint == 2, ]
  
  #helpful info for stats (n per group)
  nUnit = c(sum(testData$group == 1), sum(testData$group ==2)) 
  
  ## loop over behavior types 
for (bt in behavTypes) {
    ## set up helpful info for the model
    descript = bt
    comparison = "Flicker condition: 40 Hz vs. Random"
    testunit = "trials"
    significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
    
    ## model, anova
    obj.lmer1 = lmerTest::lmer(get(bt) ~ group + (1|animal/day), data = testData)
    lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
    
    ## set up strings for the output
    `Group Size` = paste("40 Hz: n =", as.character(nUnit[1]), "; Random: n =", as.character(nUnit[2]))
    analysisUnit = "trials"
    `P-value` = round(lmer1anova$`Pr(>F)`, digits = 5)
    `F-value` = paste("F", as.character(round(lmer1anova$NumDF, digits = 5)), ",", as.character(round(lmer1anova$DenDF, digits = 5)), " = ", as.character(round(lmer1anova$`F value`, digits = 5)), sep="")
    
    ## compile into output frame
    tmpStats = data.frame(descript, analysisUnit, comparison, `Group Size`, `P-value`, significance,`F-value`)
    outputStats = rbind(outputStats, tmpStats)
    
    ## output all stats details
    statStr = paste(descript, "|", comparison, "|", "n =", as.character(nUnit[1]), ", n =", as.character(nUnit[2]), testunit, "|", "p =", as.character(round(lmer1anova$`Pr(>F)`, digits = 5)), "|", "F", as.character(round(lmer1anova$NumDF, digits = 5)), as.character(round(lmer1anova$DenDF, digits = 5)), "=", as.character(round(lmer1anova$`F value`, digits = 5)))
    
    sink(txtFilename, append = TRUE)
    print(statStr)
    print(summary(obj.lmer1))
    print(anova(obj.lmer1, ddf="Kenward-Roger"))
    print("------------------------------------------------------------------------")
    sink()
    closeAllConnections()
    
}
  
write.csv(outputStats, csvFilename)
  
  ##### BEHAVIOR - DAY ####
  ## load data file
  fileN = "TableData_Behavior_Day.txt"
  dt = "Behavior_Day"
  behavTypes <-c("mean_duration", "mean_lickDI", "mean_nLicks", "mean_licklatency_s", "mean_licklatency_pos", "mean_trial_speed", "mean_AZ_speed", "mean_RZ_speed", "mean_Ctrl_speed", "mean_AZ_speed_slope", "mean_AZ_lick_slope")
  
  
  outputStats = data.frame()
  
  tblFilename = paste0(foldN, fileN)
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
  csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table
  
  # initialize the file
  close(file(txtFilename, open="w" ) )
  
  #read the data
  mydata <- read.csv(tblFilename, head=T)
  
  #factorize variables to be treated as categories/levels instead of numerical value
  mydata$day = as.factor(mydata$day)
  mydata$mean_animal = as.factor(mydata$mean_animal)
  mydata$mean_group = as.factor(mydata$mean_group)
  
  testData = mydata[mydata$mean_timepoint == 2,]
  
  #helpful info for stats (n per group)
  nUnit = c(sum(testData$mean_group == 1), sum(testData$mean_group ==2)) 
  
  ## loop over all behavior stats
  for (bt in behavTypes) {
    ## set up helpful info for the model
    descript = bt
    comparison = "Flicker condition: 40 Hz vs. Random"
    testunit = "days"
    significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file
    
    ## model, anova
    obj.lmer1 = lmerTest::lmer(get(bt) ~ mean_group + (1|mean_animal), data = testData)
    lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
    
    ## set up strings for the output
    `Group Size` = paste("40 Hz: n =", as.character(nUnit[1]), "; Random: n =", as.character(nUnit[2]))
    analysisUnit = "trials"
    `P-value` = round(lmer1anova$`Pr(>F)`, digits = 5)
    `F-value` = paste("F", as.character(round(lmer1anova$NumDF, digits = 5)), ",", as.character(round(lmer1anova$DenDF, digits = 5)), " = ", as.character(round(lmer1anova$`F value`, digits = 5)), sep="")
    
    ## compile into output frame
    tmpStats = data.frame(descript, analysisUnit, comparison, `Group Size`, `P-value`, significance,`F-value`)
    outputStats = rbind(outputStats, tmpStats)
    
    ## output all stats details
    statStr = paste(descript, "|", comparison, "|", "n =", as.character(nUnit[1]), ", n =", as.character(nUnit[2]), testunit, "|", "p =", as.character(round(lmer1anova$`Pr(>F)`, digits = 5)), "|", "F", as.character(round(lmer1anova$NumDF, digits = 5)), as.character(round(lmer1anova$DenDF, digits = 5)), "=", as.character(round(lmer1anova$`F value`, digits = 5)))
    
    sink(txtFilename, append = TRUE)
    print(statStr)
    print(summary(obj.lmer1))
    print(anova(obj.lmer1, ddf="Kenward-Roger"))
    print("------------------------------------------------------------------------")
    sink()
    closeAllConnections()
  }
  
  write.csv(outputStats, csvFilename)
  
  ##### PLACE CELL PROPERTIES - D2R ##### 
  ## load data file
  fileN = "TableData_PlaceCellProperties_d2r.txt"
  dt = "placecellproperties"
  propTypes <-c("SI", "sparsity", "peakpos", "maxFR", "meanFR")
  
  #initialize table structure to add stats to
  outputStats = data.frame()
  
  tblFilename = paste0(foldN, fileN)
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
  csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table
  
  #initialize the text file
  close(file(txtFilename, open="w" ) )
  
  #read data
  mydata <- read.csv(tblFilename, head=T)
  
  #factorize variables to be treated as categories/levels instead of numerical value
  mydata$day = as.factor(mydata$day)
  mydata$animal = as.factor(mydata$animal)
  mydata$group = as.factor(mydata$group)
  
  ## get data subset of interest
  testData <- mydata[mydata$isPC == 1 & mydata$timepoint == 'post',] # only place cells
  
  #helpful info for stats (n per group)
  nUnit = c(sum(testData$group == 'gamma'), sum(testData$group == 'random')) 
  
  ## loop over property types 
  for (pt in propTypes) {
    print(pt)
    ## set up helpful info for the model
    descript = pt
    comparison = "Flicker condition: 40 Hz vs. Random"
    testunit = "cells"
    significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
    
    ## model, anova
    obj.lmer1 = lmerTest::lmer(get(pt) ~ group + (1|animal), data = testData)
    lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
    
    ## set up strings for the output
    `Group Size` = paste("40 Hz: n =", as.character(nUnit[1]), "; Random: n =", as.character(nUnit[2]))
    analysisUnit = "trials"
    `P-value` = round(lmer1anova$`Pr(>F)`, digits = 5)
    `F-value` = paste("F", as.character(round(lmer1anova$NumDF, digits = 5)), ",", as.character(round(lmer1anova$DenDF, digits = 5)), " = ", as.character(round(lmer1anova$`F value`, digits = 5)), sep="")
    
    ## compile into output frame
    tmpStats = data.frame(descript, testunit, comparison, `Group Size`, `P-value`, significance,`F-value`)
    outputStats = rbind(outputStats, tmpStats)
    
    ## output all stats details
    statStr = paste(descript, "|", comparison, "|", "n =", as.character(nUnit[1]), ", n =", as.character(nUnit[2]), testunit, "|", "p =", as.character(round(lmer1anova$`Pr(>F)`, digits = 5)), "|", "F", as.character(round(lmer1anova$NumDF, digits = 5)), as.character(round(lmer1anova$DenDF, digits = 5)), "=", as.character(round(lmer1anova$`F value`, digits = 5)))
    
    sink(txtFilename, append = TRUE)
    print(statStr)
    print(summary(obj.lmer1))
    print(anova(obj.lmer1, ddf="Kenward-Roger"))
    print("-----------------------------------------------------------")
    sink()
    closeAllConnections()
  }
  
  write.csv(outputStats, csvFilename)
  
  ##### RIPPLE PROPERTIES ##### 
  #ALP 10/10/23 this is for all rippples in any location
  ## load data file
  fileN = "TableData_RippleProperties_SWR.txt"
  dt = "rippleproperties_swr"
  propTypes <-c("rippleSize", "rippleDurS")
  
  #initialize table structure to add stats to
  outputStats = data.frame()
  
  tblFilename = paste0(foldN, fileN)
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
  csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table
  
  # initialize text file
  close(file(txtFilename, open="w" ) )
  
  #read data
  mydata <- read.csv(tblFilename, head=T)
  
  #factorize variables to be treated as categories/levels instead of numerical value
  mydata$day = as.factor(mydata$day)
  mydata$animal = as.factor(mydata$animal)
  mydata$group = as.factor(mydata$group)
  
  ## get data subset of interest
  testData <- mydata[mydata$timepoint == 'post',] 
  
  #helpful info for stats (n per group)
  nUnit = c(sum(testData$group == 'gamma'), sum(testData$group == 'random')) 
  
  ## loop over property types 
  for (pt in propTypes) {
    print(pt)
    ## set up helpful info for the model
    descript = pt
    comparison = "Flicker condition: 40 Hz vs. Random"
    testunit = "ripples"
    significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
    
    ## model, anova
    obj.lmer1 = lmerTest::lmer(get(pt) ~ group + (1|day), data = testData)
    lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
    
    ## set up strings for the output
    `Group Size` = paste("40 Hz: n =", as.character(nUnit[1]), "; Random: n =", as.character(nUnit[2]))
    analysisUnit = "trials"
    `P-value` = round(lmer1anova$`Pr(>F)`, digits = 5)
    `F-value` = paste("F", as.character(round(lmer1anova$NumDF, digits = 5)), ",", as.character(round(lmer1anova$DenDF, digits = 5)), " = ", as.character(round(lmer1anova$`F value`, digits = 5)), sep="")
    
    ## compile into output frame
    tmpStats = data.frame(descript, testunit, comparison, `Group Size`, `P-value`, significance,`F-value`)
    outputStats = rbind(outputStats, tmpStats)
    
    ## output all stats details
    statStr = paste(descript, "|", comparison, "|", "n =", as.character(nUnit[1]), ", n =", as.character(nUnit[2]), testunit, "|", "p =", as.character(round(lmer1anova$`Pr(>F)`, digits = 5)), "|", "F", as.character(round(lmer1anova$NumDF, digits = 5)), as.character(round(lmer1anova$DenDF, digits = 5)), "=", as.character(round(lmer1anova$`F value`, digits = 5)))
    
    sink(txtFilename, append = TRUE)
    print(statStr)
    print(summary(obj.lmer1))
    print(anova(obj.lmer1, ddf="Kenward-Roger"))
    print("-------------------------------------------------------------------")
    sink()
    closeAllConnections()
  }
  
  write.csv(outputStats, csvFilename)
  
  
  ##### RIPPLE PROPERTIES - duration ##### 
  #ALP 10/10/23 this is for all ripples in any position
  ## load data file
  fileN = "TableData_RippleProperties_SWR.txt"
  dt = "rippleproperties_swr20NT"
  propTypes <-c("rippleDurS")
  
  #initialize table structure to add stats to
  outputStats = data.frame()
  
  tblFilename = paste0(foldN, fileN)
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
  csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table
  
  # initialize text file
  close(file(txtFilename, open="w" ) )
  
  #read data
  mydata <- read.csv(tblFilename, head=T)
  
  ## get data subset of interest
  testData <- mydata[mydata$timepoint == 'post' & mydata$perDurByRip > 20,] 
  #testData <- mydata
  
  #factorize variables to be treated as categories/levels instead of numerical value
  testData$day = as.factor(testData$day)
  testData$animal = as.factor(testData$animal)
  testData$group = as.factor(testData$group)
  testData$timepoint = as.factor(testData$timepoint)
  
  #helpful info for stats (n per group)
  nUnit = c(sum(testData$group == 'gamma'), sum(testData$group == 'random')) 
  
  ## loop over property types 
    ## set up helpful info for the model
    descript = 'ripple duration'
    comparison = "Flicker condition: 40 Hz vs. Random"
    testunit = "ripples"
    significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
    
    ## model, anova
    obj.lmer1 = lmerTest::lmer(rippleDurS ~ group + (1|animal), data = testData)
    lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
    
    ## set up strings for the output
    `Group Size` = paste("40 Hz: n =", as.character(nUnit[1]), "; Random: n =", as.character(nUnit[2]))
    analysisUnit = "trials"
    `P-value` = round(lmer1anova$`Pr(>F)`, digits = 5)
    `F-value` = paste("F", as.character(round(lmer1anova$NumDF, digits = 5)), ",", as.character(round(lmer1anova$DenDF, digits = 5)), " = ", as.character(round(lmer1anova$`F value`, digits = 5)), sep="")
    
    ## compile into output frame
    tmpStats = data.frame(descript, testunit, comparison, `Group Size`, `P-value`, significance,`F-value`)
    outputStats = rbind(outputStats, tmpStats)
    
    ## output all stats details
    statStr = paste(descript, "|", comparison, "|", "n =", as.character(nUnit[1]), ", n =", as.character(nUnit[2]), testunit, "|", "p =", as.character(round(lmer1anova$`Pr(>F)`, digits = 5)), "|", "F", as.character(round(lmer1anova$NumDF, digits = 5)), as.character(round(lmer1anova$DenDF, digits = 5)), "=", as.character(round(lmer1anova$`F value`, digits = 5)))
    
    sink(txtFilename)
    print("ripple duration - NT periods > 20s")
    print(report(obj.lmer1))
    print(statStr)
    print(summary(obj.lmer1))
    print(anova(obj.lmer1, ddf="Kenward-Roger"))
    sink()
    closeAllConnections()
  
  write.csv(outputStats, csvFilename)
  
  ##### NONTHETA PROPERTIES - ripple abundance ##### 
  # ALP 10/10/23 this is for all ripples in any position
  ## load data file
  fileN = "TableData_RippleProperties_NT.txt"
  dt = "rippleproperties_abundance"
  propTypes <-c("rippleRatePeriod")
  
  #initialize table structure to add stats to
  outputStats = data.frame()
  
  tblFilename = paste0(foldN, fileN)
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
  csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table
  
  # initialize text file
  close(file(txtFilename, open="w" ) )
  
  #read data
  mydata <- read.csv(tblFilename, head=T)
  
  ## get data subset of interest
  testData <- mydata[mydata$timepoint == 'post' & mydata$nonThetaDurS > 5, ] 

  #factorize variables to be treated as categories/levels instead of numerical value
  testData$day = as.factor(testData$day)
  testData$animal = as.factor(testData$animal)
  testData$group = as.factor(testData$group)
  testData$timepoint = as.factor(testData$timepoint)
  
  #helpful info for stats (n per group)
  nUnit = c(sum(testData$group == 'gamma'), sum(testData$group == 'random')) 
  
  ## loop over property types 
  ## set up helpful info for the model
  descript = 'ripple duration'
  comparison = "Flicker condition: 40 Hz vs. Random"
  testunit = "ripples"
  significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
  
  ## model, anova
  obj.lmer1 = lmerTest::lmer(rippleRatePeriod ~ group + (1|animal), data = testData)
  lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
  
  ## set up strings for the output
  `Group Size` = paste("40 Hz: n =", as.character(nUnit[1]), "; Random: n =", as.character(nUnit[2]))
  analysisUnit = "trials"
  `P-value` = round(lmer1anova$`Pr(>F)`, digits = 5)
  `F-value` = paste("F", as.character(round(lmer1anova$NumDF, digits = 5)), ",", as.character(round(lmer1anova$DenDF, digits = 5)), " = ", as.character(round(lmer1anova$`F value`, digits = 5)), sep="")
  
  ## compile into output frame
  tmpStats = data.frame(descript, testunit, comparison, `Group Size`, `P-value`, significance,`F-value`)
  outputStats = rbind(outputStats, tmpStats)
  
  ## output all stats details
  statStr = paste(descript, "|", comparison, "|", "n =", as.character(nUnit[1]), ", n =", as.character(nUnit[2]), testunit, "|", "p =", as.character(round(lmer1anova$`Pr(>F)`, digits = 5)), "|", "F", as.character(round(lmer1anova$NumDF, digits = 5)), as.character(round(lmer1anova$DenDF, digits = 5)), "=", as.character(round(lmer1anova$`F value`, digits = 5)))
  
  sink(txtFilename)
  print("ripple duration - NT periods > 20s")
  print(report(obj.lmer1))
  print(statStr)
  print(summary(obj.lmer1))
  print(anova(obj.lmer1, ddf="Kenward-Roger"))
  sink()
  closeAllConnections()
  
  write.csv(outputStats, csvFilename)
  
  
  ##### RIPPLE PROPERTIES - Per Day##### 
  # ALP 10/10/23 this is ripple properties for all ripples across the track
  ## load data file
  fileN = "TableData_RippleProperties_Day.txt"
  dt = "rippleproperties_day"
  propTypes <-c("mean_rippleSize", "mean_rippleDurS")
  
  #initialize table structure to add stats to
  outputStats = data.frame()
  
  tblFilename = paste0(foldN, fileN)
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
  csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table
  
  # initialize table
  close(file(txtFilename, open="w" ) )
  
  #read data
  mydata <- read.csv(tblFilename, head=T)
  
  #factorize variables to be treated as categories/levels instead of numerical value
  mydata$day = as.factor(mydata$day)
  mydata$animal = as.factor(mydata$animal)
  mydata$group = as.factor(mydata$group)
  
  ## get data subset of interest
  testData <- mydata[mydata$timepoint == 'post',] 
  
  #helpful info for stats (n per group)
  nUnit = c(sum(testData$group == 'gamma'), sum(testData$group == 'random')) 
  
  ## loop over property types 
  for (pt in propTypes) {
    print(pt)
    ## set up helpful info for the model
    descript = pt
    comparison = "Flicker condition: 40 Hz vs. Random"
    testunit = "day"
    significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
    
    ## model, anova
    obj.lmer1 = lmerTest::lmer(get(pt) ~ group + (1|animal), data = testData)
    lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
    
    ## set up strings for the output
    `Group Size` = paste("40 Hz: n =", as.character(nUnit[1]), "; Random: n =", as.character(nUnit[2]))
    analysisUnit = "trials"
    `P-value` = round(lmer1anova$`Pr(>F)`, digits = 5)
    `F-value` = paste("F", as.character(round(lmer1anova$NumDF, digits = 5)), ",", as.character(round(lmer1anova$DenDF, digits = 5)), " = ", as.character(round(lmer1anova$`F value`, digits = 5)), sep="")
    
    ## compile into output frame
    tmpStats = data.frame(descript, testunit, comparison, `Group Size`, `P-value`, significance,`F-value`)
    outputStats = rbind(outputStats, tmpStats)
    
    ## output all stats details
    statStr = paste(descript, "|", comparison, "|", "n =", as.character(nUnit[1]), ", n =", as.character(nUnit[2]), testunit, "|", "p =", as.character(round(lmer1anova$`Pr(>F)`, digits = 5)), "|", "F", as.character(round(lmer1anova$NumDF, digits = 5)), as.character(round(lmer1anova$DenDF, digits = 5)), "=", as.character(round(lmer1anova$`F value`, digits = 5)))
    
    sink(txtFilename, append = TRUE)
    print(statStr)
    print(summary(obj.lmer1))
    print(anova(obj.lmer1, ddf="Kenward-Roger"))
    print("--------------------------------------------")
    sink()
    closeAllConnections()
  }
  
  write.csv(outputStats, csvFilename)
  
  
  ##### THETA SEQUENCES - PER TRIAL ##### 
  
  ## load data file
  fileN = "TableData_ThetaSequences_PerTrial.txt"
  dt = "thetaseq_pertrial"
  propTypes <-c("PCR_loc_trial")
  
  #initialize table structure to add stats to
  outputStats = data.frame()
  
  tblFilename = paste0(foldN, fileN)
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
  csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table
  
  #initialize text file
  close(file(txtFilename, open="w" ) )
  
  #read data
  mydata <- read.csv(tblFilename, head=T)
  
  ## get data subset of interest
  testData <- mydata[mydata$timepoint == 'post' & mydata$engaged == 1 & mydata$fullTrial == 1 & mydata$rewarded == 1 & is.nan(mydata$PCR_loc_trial) == FALSE,] 
  #testData <- mydata
  
  #factorize variables to be treated as categories/levels instead of numerical value
  testData$day = as.factor(testData$day)
  testData$animal = as.factor(testData$animal)
  testData$group = as.factor(testData$group)
  testData$timepoint = as.factor(testData$timepoint)
  
  #helpful info for stats (n per group)
  nUnit = c(sum(testData$group == 'gamma'), sum(testData$group == 'random')) 
  
  ## loop over property types 
  ## set up helpful info for the model
  descript = 'RRZ theta prospective coding per trial'
  comparison = "Flicker condition: 40 Hz vs. Random"
  testunit = "trials"
  significance = '=IF(F2<0.0001,"****",IF(F2<0.001,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
  
  ## model, anova
  obj.lmer1 = lmerTest::lmer(PCR_loc_trial ~ group + (1|animal), data = testData)
  lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
  
  ## set up strings for the output
  `Group Size` = paste("40 Hz: n =", as.character(nUnit[1]), "; Random: n =", as.character(nUnit[2]))
  analysisUnit = "trials"
  `P-value` = round(lmer1anova$`Pr(>F)`, digits = 5)
  `F-value` = paste("F", as.character(round(lmer1anova$NumDF, digits=5)), ",", as.character(round(lmer1anova$DenDF, digits = 5)), " = ", as.character(round(lmer1anova$`F value`, digits = 5)), sep="")
  
  ## compile into output frame
  tmpStats = data.frame(descript, testunit, comparison, `Group Size`, `P-value`, significance,`F-value`)
  outputStats = rbind(outputStats, tmpStats)
  
  ## output all stats details
  statStr = paste(descript, "|", comparison, "|", "n =", as.character(nUnit[1]), ", n =", as.character(nUnit[2]), testunit, "|", "p =", as.character(round(lmer1anova$`Pr(>F)`, digits = 5)), "|", "F", as.character(round(lmer1anova$NumDF, digits = 5)), as.character(round(lmer1anova$DenDF, digits = 5)), "=", as.character(round(lmer1anova$`F value`, digits = 5)))
  
  sink(txtFilename, append = TRUE)
  print(report(obj.lmer1))
  print(statStr)
  print(summary(obj.lmer1))
  print(anova(obj.lmer1, ddf="Kenward-Roger"))
  sink()
  closeAllConnections()
  
  write.csv(outputStats, csvFilename)
  
  
  ##### THETA SEQUENCES - PER DAY ##### 
  
  ## load data file
  fileN = "TableData_ThetaSequences_PerDay.txt"
  dt = "thetaseq_perday"
  propTypes <-c("mean_PCR_loc_trial")
  
  #initialize table structure to add stats to
  outputStats = data.frame()
  
  tblFilename = paste0(foldN, fileN)
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
  csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table
  
  # initialize text file
  close(file(txtFilename, open="w" ) )
  
  #read data
  mydata <- read.csv(tblFilename, head=T)
  
  ## get data subset of interest
  testData <- mydata[mydata$timepoint == 'post',] 

  #factorize variables to be treated as categories/levels instead of numerical value
  testData$day = as.factor(testData$day)
  testData$animal = as.factor(testData$animal)
  testData$group = as.factor(testData$group)
  testData$timepoint = as.factor(testData$timepoint)
  
  #helpful info for stats (n per group)
  nUnit = c(sum(testData$group == 'gamma'), sum(testData$group == 'random')) 
  
  ## loop over property types 
  ## set up helpful info for the model
  descript = 'RRZ theta prospective coding per trial'
  comparison = "Flicker condition: 40 Hz vs. Random"
  testunit = "day"
  significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
  
  ## model, anova
  obj.lmer1 = lmerTest::lmer(mean_PCR_loc_trial ~ group + (1|animal), data = testData)
  lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
  
  ## set up strings for the output
  `Group Size` = paste("40 Hz: n =", as.character(nUnit[1]), "; Random: n =", as.character(nUnit[2]))
  analysisUnit = "trials"
  `P-value` = round(lmer1anova$`Pr(>F)`, digits = 5)
  `F-value` = paste("F", as.character(round(lmer1anova$NumDF, digits=5)), ",", as.character(round(lmer1anova$DenDF, digits = 5)), " = ", as.character(round(lmer1anova$`F value`, digits = 5)), sep="")
  
  ## compile into output frame
  tmpStats = data.frame(descript, testunit, comparison, `Group Size`, `P-value`, significance,`F-value`)
  outputStats = rbind(outputStats, tmpStats)
  
  ## output all stats details
  statStr = paste(descript, "|", comparison, "|", "n =", as.character(nUnit[1]), ", n =", as.character(nUnit[2]), testunit, "|", "p =", as.character(round(lmer1anova$`Pr(>F)`, digits = 5)), "|", "F", as.character(round(lmer1anova$NumDF, digits=5)), as.character(round(lmer1anova$DenDF, digits = 5)), "=", as.character(round(lmer1anova$`F value`, digits = 51)))
  
  sink(txtFilename, append = TRUE)
  print(report(obj.lmer1))
  print(statStr)
  print(summary(obj.lmer1))
  print(anova(obj.lmer1, ddf="Kenward-Roger"))
  sink()
  closeAllConnections()
  
  write.csv(outputStats, csvFilename)
  ##### THETA SEQUENCES - SEQUENCE STRENGTH - PER DAY ##### 
  
  ## load data file
  fileN = "TableData_ThetaSequences_PerDay.txt"
  dt = "thetaseq_seqstrength_perday"
  propTypes <-c("mean_dayQRatio")
  
  #initialize text file
  close(file(txtFilename, open="w" ) )
  
  #initialize table structure to add stats to
  outputStats = data.frame()
  
  tblFilename = paste0(foldN, fileN)
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
  csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table
  
  #read data
  mydata <- read.csv(tblFilename, head=T)
  
  
  ## get data subset of interest
  testData <- mydata[mydata$timepoint == 'post',] 
  
  #factorize variables to be treated as categories/levels instead of numerical value
  testData$day = as.factor(testData$day)
  testData$animal = as.factor(testData$animal)
  testData$group = as.factor(testData$group)
  testData$timepoint = as.factor(testData$timepoint)
  
  #helpful info for stats (n per group)
  nUnit = c(sum(testData$group == 'gamma'), sum(testData$group == 'random')) 
  
  ## loop over property types 
  ## set up helpful info for the model
  descript = 'RRZ theta prospective coding per trial'
  comparison = "Flicker condition: 40 Hz vs. Random"
  testunit = "day"
  significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
  
  ## model, anova
  obj.lmer1 = lmerTest::lmer(mean_dayQRatio ~ group + (1|animal), data = testData)
  lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
  
  ## set up strings for the output
  `Group Size` = paste("40 Hz: n =", as.character(nUnit[1]), "; Random: n =", as.character(nUnit[2]))
  analysisUnit = "trials"
  `P-value` = round(lmer1anova$`Pr(>F)`, digits = 5)
  `F-value` = paste("F", as.character(round(lmer1anova$NumDF, digits=5)), ",", as.character(round(lmer1anova$DenDF, digits = 5)), " = ", as.character(round(lmer1anova$`F value`, digits = 5)), sep="")
  
  ## compile into output frame
  tmpStats = data.frame(descript, testunit, comparison, `Group Size`, `P-value`, significance,`F-value`)
  outputStats = rbind(outputStats, tmpStats)
  
  ## output all stats details
  statStr = paste(descript, "|", comparison, "|", "n =", as.character(nUnit[1]), ", n =", as.character(nUnit[2]), testunit, "|", "p =", as.character(round(lmer1anova$`Pr(>F)`, digits = 5)), "|", "F", as.character(round(lmer1anova$NumDF, digits=5)), as.character(round(lmer1anova$DenDF, digits = 5)), "=", as.character(round(lmer1anova$`F value`, digits = 51)))
  
  sink(txtFilename, append = TRUE)
  print(report(obj.lmer1))
  print(statStr)
  print(summary(obj.lmer1))
  print(anova(obj.lmer1, ddf="Kenward-Roger"))
  sink()
  closeAllConnections()
  
  write.csv(outputStats, csvFilename)
  
  ##### Ripple Decoding - PER Ripple ##### 
  
  ## load data file
  fileN = "TableData_RippleDecoding_SWR.txt"
  dt = "rippledecoding_rewardzone_perSWR"
  propTypes <-c("rippleRatio", "rippleRelPos")
  
  #initialize table structure to add stats to
  outputStats = data.frame()
  tblFilename = paste0(foldN, fileN)
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
  csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table
  
  #initialize text file
  close(file(txtFilename, open="w" ) )
  
  #read data
  mydata <- read.csv(tblFilename, head=T)
  
  ## get data subset of interest
  testData <- mydata[mydata$timepoint == 'post' & mydata$includeRipple == 1 & mydata$includeDay == 1 & mydata$RZRipple == 1,] 
  #testData <- mydata
  
  #factorize variables to be treated as categories/levels instead of numerical value
  testData$day = as.factor(testData$day)
  testData$animal = as.factor(testData$animal)
  testData$group = as.factor(testData$group)
  testData$timepoint = as.factor(testData$timepoint)
  
  #helpful info for stats (n per group)
  nUnit = c(sum(testData$group == 'gamma'), sum(testData$group == 'random')) 
  
  ## loop over property types 
  for (pt in propTypes) {
  ## set up helpful info for the model
    txtFilename = paste0(foldN, "/stats_lme4_", dt, "_", pt, ".txt") #create text file to save LMM outputs
  descript = pt
  comparison = "Flicker condition: 40 Hz vs. Random"
  testunit = "ripple"
  significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
  
  ## model, anova
  obj.lmer1 = lmerTest::lmer(get(pt) ~ group + (1|animal), data = testData)
  lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
  
  ## set up strings for the output
  `Group Size` = paste("40 Hz: n =", as.character(nUnit[1]), "; Random: n =", as.character(nUnit[2]))
  analysisUnit = "trials"
  `P-value` = round(lmer1anova$`Pr(>F)`, digits = 4)
  `F-value` = paste("F", as.character(round(lmer1anova$NumDF, digits=5)), ",", as.character(round(lmer1anova$DenDF, digits = 5)), " = ", as.character(round(lmer1anova$`F value`, digits = 5)), sep="")
  
  ## compile into output frame
  tmpStats = data.frame(descript, testunit, comparison, `Group Size`, `P-value`, significance,`F-value`)
  outputStats = rbind(outputStats, tmpStats)
  
  ## output all stats details
  statStr = paste(descript, "|", comparison, "|", "n =", as.character(nUnit[1]), ", n =", as.character(nUnit[2]), testunit, "|", "p =", as.character(round(lmer1anova$`Pr(>F)`, digits = 5)), "|", "F", as.character(round(lmer1anova$NumDF, digits = 5)), as.character(round(lmer1anova$DenDF, digits = 5)), "=", as.character(round(lmer1anova$`F value`, digits = 5)))
  
  sink(txtFilename)
  print(statStr)
  print(summary(obj.lmer1))
  print(anova(obj.lmer1, ddf="Kenward-Roger"))
  sink()
  closeAllConnections()
  }
  write.csv(outputStats, csvFilename)
  
  ##### COACTIVATION DURING RIPPLES - PER PC PAIR ##### 
  #ALP 10/10/23 this is for all positions
  ## load data file
  fileN = "TableData_ripples_coactivation_placecells.txt"
  dt = "ripplecoactivation_placecells_pairs"
  propTypes <-c("fracCoactive")
  
  #initialize table structure to add stats to
  outputStats = data.frame()
  
  tblFilename = paste0(foldN, fileN)
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
  csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table
  
  #read data
  mydata <- read.csv(tblFilename, head=T)
  
  ## get data subset of interest
  testData <- mydata[mydata$timepoint == 'post' & mydata$nDayRipples > 10,] 

  #factorize variables to be treated as categories/levels instead of numerical value
  testData$day = as.factor(testData$day)
  testData$animal = as.factor(testData$animal)
  testData$group = as.factor(testData$group)
  testData$timepoint = as.factor(testData$timepoint)
  
  #helpful info for stats (n per group)
  nUnit = c(sum(testData$group == 'gamma'), sum(testData$group == 'random')) 
  
  ## loop over property types 
  for (pt in propTypes) {
    ## set up helpful info for the model
    descript = pt
    comparison = "Flicker condition: 40 Hz vs. Random"
    testunit = "place cell pairs"
    significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
    
    ## model, anova
    obj.lmer1 = lmerTest::lmer(get(pt) ~ group + (1|animal), data = testData)
    lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
    
    ## set up strings for the output
    `Group Size` = paste("40 Hz: n =", as.character(nUnit[1]), "; Random: n =", as.character(nUnit[2]))
    analysisUnit = "trials"
    `P-value` = round(lmer1anova$`Pr(>F)`, digits = 5)
    `F-value` = paste("F", as.character(round(lmer1anova$NumDF, digits = 5)), ",", as.character(round(lmer1anova$DenDF, digits = 5)), " = ", as.character(round(lmer1anova$`F value`, digits = 5)), sep="")
    
    ## compile into output frame
    tmpStats = data.frame(descript, testunit, comparison, `Group Size`, `P-value`, significance,`F-value`)
    outputStats = rbind(outputStats, tmpStats)
    
    ## output all stats details
    statStr = paste(descript, "|", comparison, "|", "n =", as.character(nUnit[1]), ", n =", as.character(nUnit[2]), testunit, "|", "p =", as.character(round(lmer1anova$`Pr(>F)`, digits = 5)), "|", "F", as.character(round(lmer1anova$NumDF, digits = 5)), as.character(round(lmer1anova$DenDF, digits = 5)), "=", as.character(round(lmer1anova$`F value`, digits = 5)))
    
    sink(txtFilename)
    print(statStr)
    print(summary(obj.lmer1))
    print(anova(obj.lmer1, ddf="Kenward-Roger"))
    sink()
    closeAllConnections()
  }
  write.csv(outputStats, csvFilename)
  
  ##### ACTIVATION DURING RIPPLES - PER PC PAIR ##### 
  #ALP 10/10/23 this is for all positions
  ## load data file
  fileN = "TableData_ripples_activation_placecells.txt"
  dt = "rippleactivation_placecells_pairs"
  propTypes <-c("fracActive")
  
  #initialize table structure to add stats to
  outputStats = data.frame()
  
  tblFilename = paste0(foldN, fileN)
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
  csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table
  
  #read data
  mydata <- read.csv(tblFilename, head=T)
  
  ## get data subset of interest
  testData <- mydata[mydata$timepoint == 'post' & mydata$nDayRipples > 10, ] 
  
  #factorize variables to be treated as categories/levels instead of numerical value
  testData$day = as.factor(testData$day)
  testData$animal = as.factor(testData$animal)
  testData$group = as.factor(testData$group)
  testData$timepoint = as.factor(testData$timepoint)
  
  #helpful info for stats (n per group)
  nUnit = c(sum(testData$group == 'gamma'), sum(testData$group == 'random')) 
  
  ## loop over property types 
  for (pt in propTypes) {
    ## set up helpful info for the model
    descript = pt
    comparison = "Flicker condition: 40 Hz vs. Random"
    testunit = "place cell pairs"
    significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
    
    ## model, anova
    obj.lmer1 = lmerTest::lmer(get(pt) ~ group + (1|animal), data = testData)
    lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
    
    ## set up strings for the output
    `Group Size` = paste("40 Hz: n =", as.character(nUnit[1]), "; Random: n =", as.character(nUnit[2]))
    analysisUnit = "trials"
    `P-value` = round(lmer1anova$`Pr(>F)`, digits = 5)
    `F-value` = paste("F", as.character(round(lmer1anova$NumDF, digits = 5)), ",", as.character(round(lmer1anova$DenDF, digits = 5)), " = ", as.character(round(lmer1anova$`F value`, digits = 5)), sep="")
    
    ## compile into output frame
    tmpStats = data.frame(descript, testunit, comparison, `Group Size`, `P-value`, significance,`F-value`)
    outputStats = rbind(outputStats, tmpStats)
    
    ## output all stats details
    statStr = paste(descript, "|", comparison, "|", "n =", as.character(nUnit[1]), ", n =", as.character(nUnit[2]), testunit, "|", "p =", as.character(round(lmer1anova$`Pr(>F)`, digits = 5)), "|", "F", as.character(round(lmer1anova$NumDF, digits = 5)), as.character(round(lmer1anova$DenDF, digits = 5)), "=", as.character(round(lmer1anova$`F value`, digits = 5)))
    
    sink(txtFilename)
    print(statStr)
    print(summary(obj.lmer1))
    print(anova(obj.lmer1, ddf="Kenward-Roger"))
    sink()
    closeAllConnections()
  }
  write.csv(outputStats, csvFilename)
  
  ##### Ripple Decoding - Per Day ##### 
  
  ## load data file
  fileN = "TableData_RippleDecoding_SWR_Day.txt"
  dt = "rippledecoding_rewardzone_perDay"
  propTypes <-c("prospRipples","prospRipples_RZ", "nRipples", "nRipples_RZ")
  
  #initialize table structure to add stats to
  outputStats = data.frame()
  
  tblFilename = paste0(foldN, fileN)
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
  csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table
  
  #read data
  mydata <- read.csv(tblFilename, head=T)
  
  ## get data subset of interest
 testData <- mydata[mydata$timepoint == 'post' & mydata$includeDay == 1,] 

  #factorize variables to be treated as categories/levels instead of numerical value
  testData$day = as.factor(testData$day)
  testData$animal = as.factor(testData$animal)
  testData$group = as.factor(testData$group)
  testData$timepoint = as.factor(testData$timepoint)
  
  #helpful info for stats (n per group)
  nUnit = c(sum(testData$group == 'gamma'), sum(testData$group == 'random')) 
  
  ## loop over property types 
  ## set up helpful info for the model
  for (pt in propTypes) {
  descript = pt
  comparison = "Flicker condition: 40 Hz vs. Random"
  testunit = "day"
  significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
  
  ## model, anova
  obj.lmer1 = lmerTest::lmer(get(pt) ~ group + (1|animal), data = testData)
  lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
  
  ## set up strings for the output
  `Group Size` = paste("40 Hz: n =", as.character(nUnit[1]), "; Random: n =", as.character(nUnit[2]))
  analysisUnit = "trials"
  `P-value` = round(lmer1anova$`Pr(>F)`, digits = 5)
  `F-value` = paste("F", as.character(round(lmer1anova$NumDF, digits = 5)), ",", as.character(round(lmer1anova$DenDF, digits = 5)), " = ", as.character(round(lmer1anova$`F value`, digits = 5)), sep="")
  
  ## compile into output frame
  tmpStats = data.frame(descript, testunit, comparison, `Group Size`, `P-value`, significance,`F-value`)
  outputStats = rbind(outputStats, tmpStats)
  
  ## output all stats details
  statStr = paste(descript, "|", comparison, "|", "n =", as.character(nUnit[1]), ", n =", as.character(nUnit[2]), testunit, "|", "p =", as.character(round(lmer1anova$`Pr(>F)`, digits = 5)), "|", "F", as.character(round(lmer1anova$NumDF, digits = 5)), as.character(round(lmer1anova$DenDF, digits = 5)), "=", as.character(round(lmer1anova$`F value`, digits = 5)))
  
  sink(txtFilename, append = TRUE)
  #print(report(obj.lmer1))
  print(statStr)
  print(summary(obj.lmer1))
  print(anova(obj.lmer1, ddf="Kenward-Roger"))
  print ("-------------------------------------------------------")
  sink()
  closeAllConnections()
  }
  
  write.csv(outputStats, csvFilename)
  
  ##### Ripple Decoding - ALL RIPPLES - PER Ripple ##### 
  
  ## load data file
  fileN = "TableData_RippleDecoding_SWR.txt"
  dt = "rippledecoding_allzone_perSWR"
  propTypes <-c("rippleRatio", "rippleRelPos")
  
  #initialize table structure to add stats to
  outputStats = data.frame()
  tblFilename = paste0(foldN, fileN)
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
  csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table
  
  #initialize text file
  close(file(txtFilename, open="w" ) )
  
  #read data
  mydata <- read.csv(tblFilename, head=T)
  
  ## get data subset of interest
  testData <- mydata[mydata$timepoint == 'post' & mydata$includeRipple == 1 & mydata$includeDay == 1,] 
  #testData <- mydata
  
  #factorize variables to be treated as categories/levels instead of numerical value
  testData$day = as.factor(testData$day)
  testData$animal = as.factor(testData$animal)
  testData$group = as.factor(testData$group)
  testData$timepoint = as.factor(testData$timepoint)
  
  #helpful info for stats (n per group)
  nUnit = c(sum(testData$group == 'gamma'), sum(testData$group == 'random')) 
  
  ## loop over property types 
  for (pt in propTypes) {
    ## set up helpful info for the model
    descript = pt
    comparison = "Flicker condition: 40 Hz vs. Random"
    testunit = "ripple"
    significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
    
    ## model, anova
    obj.lmer1 = lmerTest::lmer(get(pt) ~ group + (1|animal), data = testData)
    lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
    
    ## set up strings for the output
    `Group Size` = paste("40 Hz: n =", as.character(nUnit[1]), "; Random: n =", as.character(nUnit[2]))
    analysisUnit = "trials"
    `P-value` = round(lmer1anova$`Pr(>F)`, digits = 4)
    `F-value` = paste("F", as.character(round(lmer1anova$NumDF, digits=5)), ",", as.character(round(lmer1anova$DenDF, digits = 5)), " = ", as.character(round(lmer1anova$`F value`, digits = 5)), sep="")
    
    ## compile into output frame
    tmpStats = data.frame(descript, testunit, comparison, `Group Size`, `P-value`, significance,`F-value`)
    outputStats = rbind(outputStats, tmpStats)
    
    ## output all stats details
    statStr = paste(descript, "|", comparison, "|", "n =", as.character(nUnit[1]), ", n =", as.character(nUnit[2]), testunit, "|", "p =", as.character(round(lmer1anova$`Pr(>F)`, digits = 5)), "|", "F", as.character(round(lmer1anova$NumDF, digits = 5)), as.character(round(lmer1anova$DenDF, digits = 5)), "=", as.character(round(lmer1anova$`F value`, digits = 5)))
    
    sink(txtFilename, append = TRUE)
    print(statStr)
    print(summary(obj.lmer1))
    print(anova(obj.lmer1, ddf="Kenward-Roger"))
    sink()
    closeAllConnections()
  }
  write.csv(outputStats, csvFilename)
  
  ##### Proportion of Goal Related PCs - Per day ##### 
  
  ## load data file
  fileN = "TableData_PlaceCellProperties_d2r_DayData.txt"
  dt = "prop_RRZPC_day"
  propTypes <-c("propRRZ_PC")
  
  #initialize table structure to add stats to
  outputStats = data.frame()
  
  tblFilename = paste0(foldN, fileN)
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
  csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table
  
  #read data
  mydata <- read.csv(tblFilename, head=T)
  
  ## get data subset of interest
  testData <- mydata[mydata$timepoint == 'post',] 
  
  #factorize variables to be treated as categories/levels instead of numerical value
  testData$day = as.factor(testData$day)
  testData$animal = as.factor(testData$animal)
  testData$group = as.factor(testData$group)
  testData$timepoint = as.factor(testData$timepoint)
  
  #helpful info for stats (n per group)
  nUnit = c(sum(testData$group == 'gamma'), sum(testData$group == 'random')) 
  
  ## loop over property types 
  ## set up helpful info for the model
  for (pt in propTypes) {
    descript = pt
    comparison = "Flicker condition: 40 Hz vs. Random"
    testunit = "day"
    significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
    
    ## model, anova
    obj.lmer1 = lmerTest::lmer(get(pt) ~ group + (1|animal), data = testData)
    lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
    
    ## set up strings for the output
    `Group Size` = paste("40 Hz: n =", as.character(nUnit[1]), "; Random: n =", as.character(nUnit[2]))
    analysisUnit = "trials"
    `P-value` = round(lmer1anova$`Pr(>F)`, digits = 5)
    `F-value` = paste("F", as.character(round(lmer1anova$NumDF, digits = 5)), ",", as.character(round(lmer1anova$DenDF, digits = 5)), " = ", as.character(round(lmer1anova$`F value`, digits = 5)), sep="")
    
    ## compile into output frame
    tmpStats = data.frame(descript, testunit, comparison, `Group Size`, `P-value`, significance,`F-value`)
    outputStats = rbind(outputStats, tmpStats)
    
    ## output all stats details
    statStr = paste(descript, "|", comparison, "|", "n =", as.character(nUnit[1]), ", n =", as.character(nUnit[2]), testunit, "|", "p =", as.character(round(lmer1anova$`Pr(>F)`, digits = 5)), "|", "F", as.character(round(lmer1anova$NumDF, digits = 5)), as.character(round(lmer1anova$DenDF, digits = 5)), "=", as.character(round(lmer1anova$`F value`, digits = 5)))
    
    sink(txtFilename, append = TRUE)
    #print(report(obj.lmer1))
    print(statStr)
    print(summary(obj.lmer1))
    print(anova(obj.lmer1, ddf="Kenward-Roger"))
    print ("-------------------------------------------------------")
    sink()
    closeAllConnections()
  }
  
  write.csv(outputStats, csvFilename)
  
  ##### THETA SEQUENCES - PCR vs. time slowed in RZ ##### 
  
  ## load data file
  fileN = "TableData_ThetaSequences_withBehavior_perTrial.txt"
  dt = "PCR_timebelowthresh"
  propTypes <-c("RRZ_time_below_thresh","time_2_nextRZ")
  
  #initialize table structure to add stats to
  outputStats = data.frame()
  
  tblFilename = paste0(foldN, fileN)
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
  csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table
  
  #read data
  mydata <- read.csv(tblFilename, head=T)
  
  ## get data subset of interest
  testData <- mydata[mydata$timepoint == 'post' & mydata$group == 'gamma'
                     & mydata$fullTrial == 1 & mydata$engaged == 1 & mydata$rewarded == 1 & mydata$significantSeq == 1
                     & mydata$PCR_bin !=0,] 
  
  #factorize variables to be treated as categories/levels instead of numerical value
  testData$day = as.factor(testData$day)
  testData$animal = as.factor(testData$animal)
  testData$group = as.factor(testData$group)
  testData$timepoint = as.factor(testData$timepoint)
  testData$PCR_bin = as.factor(testData$PCR_bin )
  
  #helpful info for stats (n per group)
  nUnit = c(sum(testData$group == 'gamma'), sum(testData$group == 'random')) 
  
  ## loop over property types 
  ## set up helpful info for the model
  for (pt in propTypes) {
    descript = pt
    comparison = "binned PCR"
    testunit = "trial"
    significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
    
    ## model, anova
    obj.lmer1 = lmerTest::lmer(get(pt) ~ PCR_bin + (1|animal), data = testData)
    lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
    
    ## set up strings for the output
    `Group Size` = paste("40 Hz: n =", as.character(nUnit[1]), "; Random: n =", as.character(nUnit[2]))
    analysisUnit = "trials"
    `P-value` = round(lmer1anova$`Pr(>F)`, digits = 5)
    `F-value` = paste("F", as.character(round(lmer1anova$NumDF, digits = 5)), ",", as.character(round(lmer1anova$DenDF, digits = 5)), " = ", as.character(round(lmer1anova$`F value`, digits = 5)), sep="")
    
    ## compile into output frame
    tmpStats = data.frame(descript, testunit, comparison, `Group Size`, `P-value`, significance,`F-value`)
    outputStats = rbind(outputStats, tmpStats)
    
    ## output all stats details
    statStr = paste(descript, "|", comparison, "|", "n =", as.character(nUnit[1]), ", n =", as.character(nUnit[2]), testunit, "|", "p =", as.character(round(lmer1anova$`Pr(>F)`, digits = 5)), "|", "F", as.character(round(lmer1anova$NumDF, digits = 5)), as.character(round(lmer1anova$DenDF, digits = 5)), "=", as.character(round(lmer1anova$`F value`, digits = 5)))
    
    sink(txtFilename, append = TRUE)
    #print(report(obj.lmer1))
    print(statStr)
    print(summary(obj.lmer1))
    print(anova(obj.lmer1, ddf="Kenward-Roger"))
    print(contrast(emmeans(obj.lmer1, specs=c("PCR_bin")), "pairwise", adjust = "tukey"))
    #print(contrast(emmeans(obj.lmer1, specs=c("PCR_bin")), method = list("bin1-bin2" = testData$PCR_bin1 - testData$PCR_bin2, "bin3-bin2" = testData$PCR_bin3 - testData$PCR_bin2), adjust = "tukey"))
    #emm1 = emmeans(obj.lmer1, specs = c("PCR_bin"))
    #print(contrast(emm1, method = list("PCR_bin1"-"PCR_bin2")))
    print ("-------------------------------------------------------")
    sink()
    closeAllConnections()
  }
  
  write.csv(outputStats, csvFilename)
  
  
  ##### THETA SEQUENCES - post reward zone PCR by control trial speed ##### 
  ## load data file
  fileN = "TableData_ProspectiveCoding_AfterReward_perTrial_postDataOnly.txt"
  dt = "thetaseq_pertrial_postRew"
  propTypes <-c("PostR_PCR_loc_trial")
  groupTypes <-c ("random", "gamma")
  
  #initialize table structure to add stats to
  outputStats = data.frame()
  
  tblFilename = paste0(foldN, fileN)
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
  csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table
  
  #initialize text file
  close(file(txtFilename, open="w" ) )
  
  #read data
  mydata <- read.csv(tblFilename, head=T)
  
  for (gt in groupTypes) {
  ## get data subset of interest
  testData <- mydata[mydata$group == gt,] 
  #testData <- mydata
  
  #factorize variables to be treated as categories/levels instead of numerical value
  
  #testData$day = as.factor(testData$day)
  testData$animal = as.factor(testData$animal)
  testData$group = as.factor(testData$group)
  #testData$timepoint = as.factor(testData$timepoint)
  testData$Ctrl_speed_subset= as.factor(testData$Ctrl_speed_subset)
  
  #helpful info for stats (n per group)
  nUnit = c(sum(testData$Ctrl_speed_subset == 1), sum(testData$Ctrl_speed_subset == 2)) 
  
  ## loop over property types 
  ## set up helpful info for the model
  descript = 'RRZ theta prospective coding per trial'
  comparison = "Speed subset: low vs. high"
  testunit = "trials"
  significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
  
  ## model, anova
  obj.lmer1 = lmerTest::lmer(PostR_PCR_loc_trial ~ Ctrl_speed_subset + (1|animal), data = testData)
  lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
  
  ## set up strings for the output
  `Group Size` = paste(gt, " -- speed 1 n =", as.character(nUnit[1]), "; speed 2 =", as.character(nUnit[2]))
  analysisUnit = "trials"
  `P-value` = round(lmer1anova$`Pr(>F)`, digits = 10)
  `F-value` = paste("F", as.character(round(lmer1anova$NumDF, digits=5)), ",", as.character(round(lmer1anova$DenDF, digits = 5)), " = ", as.character(round(lmer1anova$`F value`, digits = 5)), sep="")
  
  ## compile into output frame
  tmpStats = data.frame(descript, testunit, comparison, `Group Size`, `P-value`, significance,`F-value`)
  outputStats = rbind(outputStats, tmpStats)
  
  ## output all stats details
  statStr = paste(gt, "|", descript, "|", comparison, "|", "n =", as.character(nUnit[1]), ", n =", as.character(nUnit[2]), testunit, "|", "p =", as.character(round(lmer1anova$`Pr(>F)`, digits = 5)), "|", "F", as.character(round(lmer1anova$NumDF, digits = 5)), as.character(round(lmer1anova$DenDF, digits = 5)), "=", as.character(round(lmer1anova$`F value`, digits = 5)))
  
  sink(txtFilename, append = TRUE)
  print(paste("~~~~~~~ group = ", gt))
  print(report(obj.lmer1))
  print(statStr)
  print(summary(obj.lmer1))
  print(anova(obj.lmer1, ddf="Kenward-Roger"))
  print ("-------------------------------------------------------")
  sink()
  closeAllConnections()
  
  write.csv(outputStats, csvFilename)
  }
  
  ##### THETA SEQUENCES - MULTI COMPARE - post reward zone PCR by control trial speed ##### 
  ## load data file
  fileN = "TableData_ProspectiveCoding_AfterReward_perTrial_postDataOnly.txt"
  dt = "thetaseq_pertrial_postRew_multicompare"
  propTypes <-c("PostR_PCR_loc_trial")
  groupTypes <-c ("random", "gamma")
  
  #initialize table structure to add stats to
  outputStats = data.frame()
  
  tblFilename = paste0(foldN, fileN)
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
  csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table
  
  #initialize text file
  close(file(txtFilename, open="w" ) )
  
  #read data
  mydata <- read.csv(tblFilename, head=T)
  
  ## get data subset of interest
  #testData <- mydata[mydata$group == gt,] 
  testData <- mydata
  
  #factorize variables to be treated as categories/levels instead of numerical value
  
  #testData$day = as.factor(testData$day)
  testData$animal = as.factor(testData$animal)
  testData$group = as.factor(testData$group)
  #testData$timepoint = as.factor(testData$timepoint)
  testData$Ctrl_speed_subset= as.factor(testData$Ctrl_speed_subset)
  
  #helpful info for stats (n per group)
  nUnit = c(sum(testData$Ctrl_speed_subset == 1), sum(testData$Ctrl_speed_subset == 2)) 
  
  ## loop over property types 
  ## set up helpful info for the model
  descript = 'RRZ theta prospective coding per trial'
  comparison = "Speed subset: low vs. high"
  testunit = "trials"
  significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
  
  ## model, anova
  obj.lmer1 = lmerTest::lmer(PostR_PCR_loc_trial ~ Ctrl_speed_subset*group + (1|animal), data = testData)
  lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
  
  ## set up strings for the output
  `Group Size` = paste(gt, " -- speed 1 n =", as.character(nUnit[1]), "; speed 2 =", as.character(nUnit[2]))
  analysisUnit = "trials"
  `P-value` = round(lmer1anova$`Pr(>F)`, digits = 10)
  `F-value` = paste("F", as.character(round(lmer1anova$NumDF, digits=5)), ",", as.character(round(lmer1anova$DenDF, digits = 5)), " = ", as.character(round(lmer1anova$`F value`, digits = 5)), sep="")
  
  ## compile into output frame
  tmpStats = data.frame(descript, testunit, comparison, `Group Size`, `P-value`, significance,`F-value`)
  outputStats = rbind(outputStats, tmpStats)
  
  ## output all stats details
  statStr = paste(gt, "|", descript, "|", comparison, "|", "n =", as.character(nUnit[1]), ", n =", as.character(nUnit[2]), testunit, "|", "p =", as.character(round(lmer1anova$`Pr(>F)`, digits = 5)), "|", "F", as.character(round(lmer1anova$NumDF, digits = 5)), as.character(round(lmer1anova$DenDF, digits = 5)), "=", as.character(round(lmer1anova$`F value`, digits = 5)))
  
  sink(txtFilename, append = TRUE)
  #print(paste("~~~~~~~ group = ", gt))
  print(report(obj.lmer1))
  print(statStr)
  print(summary(obj.lmer1))
  print(anova(obj.lmer1, ddf="Kenward-Roger"))
  print(contrast(emmeans(obj.lmer1, specs=c("Ctrl_speed_subset", "group")), "pairwise", adjust = "tukey"))
  print ("-------------------------------------------------------")
  sink()
  closeAllConnections()
  
  write.csv(outputStats, csvFilename)
  
  
  
  ##### THETA SEQUENCES - RRZ PCR per trial, split by speed ##### 
  ## load data file
  fileN = "TableData_ProspectiveCoding_AfterReward_perTrial_postDataOnly.txt"
  dt = "thetaseq_pertrial_RRZPCR_speedsubsets"
  propTypes <-c("PCR_loc_trial")
  
  #initialize table structure to add stats to
  outputStats = data.frame()
  
  tblFilename = paste0(foldN, fileN)
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
  csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table
  
  #initialize text file
  close(file(txtFilename, open="w" ) )
  
  #read data
  mydata <- read.csv(tblFilename, head=T)
  for (s in 1:2){
    ## get data subset of interest
    testData <- mydata[mydata$Ctrl_speed_subset == s,] 
    
    #factorize variables to be treated as categories/levels instead of numerical value
    testData$animal = as.factor(testData$animal)
    testData$group = as.factor(testData$group)
    testData$Ctrl_speed_subset= as.factor(testData$Ctrl_speed_subset)
    
    #helpful info for stats (n per group)
    nUnit = c(sum(testData$group == 'gamma'), sum(testData$group == 'random')) 
    
    ## loop over property types 
    ## set up helpful info for the model
    descript = 'RRZ theta prospective coding per trial'
    comparison = "Flicker condition: 40 Hz vs. Random"
    testunit = "trials"
    significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
    
    ## model, anova
    obj.lmer1 = lmerTest::lmer(PCR_loc_trial ~ group + (1|animal), data = testData)
    lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
    
    ## set up strings for the output
    `Group Size` = paste("speed ", s, "; 40 Hz: n =", as.character(nUnit[1]), "; Random: n =", as.character(nUnit[2]))
    analysisUnit = "trials"
    `P-value` = round(lmer1anova$`Pr(>F)`, digits = 5)
    `F-value` = paste("F", as.character(round(lmer1anova$NumDF, digits=5)), ",", as.character(round(lmer1anova$DenDF, digits = 5)), " = ", as.character(round(lmer1anova$`F value`, digits = 5)), sep="")
    
    ## compile into output frame
    tmpStats = data.frame(descript, testunit, comparison, `Group Size`, `P-value`, significance,`F-value`)
    outputStats = rbind(outputStats, tmpStats)
    
    ## output all stats details
    statStr = paste(descript, "|", comparison, "|", "n =", as.character(nUnit[1]), ", n =", as.character(nUnit[2]), testunit, "|", "p =", as.character(round(lmer1anova$`Pr(>F)`, digits = 5)), "|", "F", as.character(round(lmer1anova$NumDF, digits = 5)), as.character(round(lmer1anova$DenDF, digits = 5)), "=", as.character(round(lmer1anova$`F value`, digits = 5)))
    
    sink(txtFilename, append = TRUE)
    print(paste("speed subset = ", s))
    print(report(obj.lmer1))
    print(statStr)
    print(summary(obj.lmer1))
    print(anova(obj.lmer1, ddf="Kenward-Roger"))
    print ("-------------------------------------------------------")
    sink()
    closeAllConnections()
    
    write.csv(outputStats, csvFilename)
  }
  
  ##### THETA SEQUENCES - control speed comparison - split by speed ##### 
  ## load data file
  fileN = "TableData_ProspectiveCoding_AfterReward_perTrial_postDataOnly.txt"
  dt = "thetaseq_pertrial_CtrlSpeed_speedsubsets"
  propTypes <-c("Ctrl_speed")
  
  #initialize table structure to add stats to
  outputStats = data.frame()
  
  tblFilename = paste0(foldN, fileN)
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
  csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table
  
  #initialize text file
  close(file(txtFilename, open="w" ) )
  
  #read data
  mydata <- read.csv(tblFilename, head=T)
  for (s in 1:2){
    ## get data subset of interest
    testData <- mydata[mydata$Ctrl_speed_subset == s,] 
    
    #factorize variables to be treated as categories/levels instead of numerical value
    testData$animal = as.factor(testData$animal)
    testData$group = as.factor(testData$group)
    testData$Ctrl_speed_subset= as.factor(testData$Ctrl_speed_subset)
    
    #helpful info for stats (n per group)
    nUnit = c(sum(testData$group == 'gamma'), sum(testData$group == 'random')) 
    
    ## loop over property types 
    ## set up helpful info for the model
    descript = 'RRZ theta prospective coding per trial'
    comparison = "Flicker condition: 40 Hz vs. Random"
    testunit = "trials"
    significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
    
    ## model, anova
    obj.lmer1 = lmerTest::lmer(Ctrl_speed ~ group + (1|animal), data = testData)
    lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
    
    ## set up strings for the output
    `Group Size` = paste("speed ", s, "; 40 Hz: n =", as.character(nUnit[1]), "; Random: n =", as.character(nUnit[2]))
    analysisUnit = "trials"
    `P-value` = round(lmer1anova$`Pr(>F)`, digits = 5)
    `F-value` = paste("F", as.character(round(lmer1anova$NumDF, digits=5)), ",", as.character(round(lmer1anova$DenDF, digits = 5)), " = ", as.character(round(lmer1anova$`F value`, digits = 5)), sep="")
    
    ## compile into output frame
    tmpStats = data.frame(descript, testunit, comparison, `Group Size`, `P-value`, significance,`F-value`)
    outputStats = rbind(outputStats, tmpStats)
    
    ## output all stats details
    statStr = paste(descript, "|", comparison, "|", "n =", as.character(nUnit[1]), ", n =", as.character(nUnit[2]), testunit, "|", "p =", as.character(round(lmer1anova$`Pr(>F)`, digits = 5)), "|", "F", as.character(round(lmer1anova$NumDF, digits = 5)), as.character(round(lmer1anova$DenDF, digits = 5)), "=", as.character(round(lmer1anova$`F value`, digits = 5)))
    
    sink(txtFilename, append = TRUE)
    print(paste("speed subset = ", s))
    print(report(obj.lmer1))
    print(statStr)
    print(summary(obj.lmer1))
    print(anova(obj.lmer1, ddf="Kenward-Roger"))
    print ("-------------------------------------------------------")
    sink()
    closeAllConnections()
    
    write.csv(outputStats, csvFilename)
  }
  
  
  ##### PPC DATA! ##### 
  ## load data file
  PPCtype <- c("VR", "Flicker")
  celltype <- c("Pyr", "Int")
  dtt <- c("ppc_crossreg_pyr", "ppc_crossreg_in")
  
  for (PPCt in PPCtype) {
    for (CT in celltype) {
      # fileN = "TableData_PPC_VR_CrossReg_Pyr.txt"
      fileN <- paste("TableData_PPC_", PPCt, "_CrossReg_", CT, ".txt", sep = "")
      dt <- paste("ppc_crossreg_", PPCt, "_", CT, sep = "")
      # propTypes <-c("VR_PYR_PPC")
      
      #initialize table structure to add stats to
      outputStats = data.frame()
      
      tblFilename = paste0(foldN, fileN)
      txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
      csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table
      
      #initialize text file
      close(file(txtFilename, open="w" ) )
      
      #read data
      mydata <- read.csv(tblFilename, head=T)
      ## get data subset of interest
      testData <- mydata[(mydata$freqID == 2 | mydata$freqID == 4 | mydata$freqID == 5) & mydata$regID == 2,] # freqID 4 = slow gamma, freqID 5 = medium gamma, regionID 2 = CA3-CA1 LFP
      
      #factorize variables to be treated as categories/levels instead of numerical value
      testData$anID = as.factor(testData$anID)
      testData$groupID = as.factor(testData$groupID)
      testData$freqID= as.factor(testData$freqID)
      testData$cellID = as.factor(testData$cellID)
      
      #helpful info for stats (n per group)
      #nUnit = c(sum(testData$group == 'gamma'), sum(testData$group == 'random')) 
      
      ## loop over property types 
      ## set up helpful info for the model
      #descript = 'RRZ theta prospective coding per trial'
      #comparison = "Flicker condition: 40 Hz vs. Random"
      #testunit = "trials"
      #significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
      
      ## model, anova
      obj.lmer1 = lmerTest::lmer(ppc ~ groupID*freqID + (1|anID/cellID), data = testData)
      lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
      
      ## set up strings for the output
      # `Group Size` = paste("speed ", s, "; 40 Hz: n =", as.character(nUnit[1]), "; Random: n =", as.character(nUnit[2]))
      # analysisUnit = "trials"
      # `P-value` = round(lmer1anova$`Pr(>F)`, digits = 5)
      # `F-value` = paste("F", as.character(round(lmer1anova$NumDF, digits=5)), ",", as.character(round(lmer1anova$DenDF, digits = 5)), " = ", as.character(round(lmer1anova$`F value`, digits = 5)), sep="")
      # 
      # ## compile into output frame
      # tmpStats = data.frame(descript, testunit, comparison, `Group Size`, `P-value`, significance,`F-value`)
      # outputStats = rbind(outputStats, tmpStats)
      # 
      # ## output all stats details
      # statStr = paste(descript, "|", comparison, "|", "n =", as.character(nUnit[1]), ", n =", as.character(nUnit[2]), testunit, "|", "p =", as.character(round(lmer1anova$`Pr(>F)`, digits = 5)), "|", "F", as.character(round(lmer1anova$NumDF, digits = 5)), as.character(round(lmer1anova$DenDF, digits = 5)), "=", as.character(round(lmer1anova$`F value`, digits = 5)))
      # 
      sink(txtFilename, append = TRUE)
      #print(paste("speed subset = ", s))
      print(report(obj.lmer1))
      #print(statStr)
      print(summary(obj.lmer1))
      print(anova(obj.lmer1, ddf="Kenward-Roger"))
      print ("-------------------------------------------------------")
      sink()
      closeAllConnections()
      
      write.csv(outputStats, csvFilename)
    }
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  