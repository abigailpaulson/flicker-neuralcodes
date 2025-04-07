### load the libraries
# ALP 10/10/23
# stats for ripple properties (reward zone only ripples) found in supplement
rm()

library(lme4)
library(emmeans)
library(sjPlot)
library(report)
library(ggplot2)

#load MATLAB output (in table format in .txt) file in R
foldN = "//ad.gatech.edu/bme/labs/singer/Abby/code/chronicflicker-ephys-prospectivecoding/results/LMM_R/"

##### RIPPLE PROPERTIES ##### 
## load data file

## load data file
fileN = "TableData_RippleDecoding_SWR.txt"
dt = "rippleproperties_swr_RZ"
propTypes <-c("rippleSize", "rippleDurationS")

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
testData <- mydata[mydata$timepoint == 'post' & mydata$RZRipple == 1,] 

#helpful info for stats (n per group)
nUnit = c(sum(testData$group == 'gamma'), sum(testData$group == 'random')) 

## loop over property types 
for (pt in propTypes) {
  print(pt)
  ## set up helpful info for the model
  descript = pt
  comparison = "Flicker condition: 40 Hz vs. Random"
  testunit = "ripples"
  significance = '=IF(F2<0.0001,"****",IF(F2<0.001,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
  
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

##### RIPPLE RATE ##### 
## load data file
fileN = "TableData_RippleRate_RZ.txt"
dt = "ripplerate_RZ"
propTypes <-c("rippleRateS")

#initialize table structure to add stats to
outputStats = data.frame()

tblFilename = paste0(foldN, fileN)
txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table

#read data
mydata <- read.csv(tblFilename, head=T)

## get data subset of interest
testData <- mydata[mydata$timepoint == 'post' & mydata$nRipples > 5,] 

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
  testunit = "days"
  significance = '=IF(F2<0.0001,"****",IF(F2<0.001,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
  
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

##### COACTIVATION DURING RIPPLES##### 
## load data file
fileN = "TableData_ripples_coactivation_placecells_RZ.txt"
dt = "ripplecoactivation_placecells_pairs_RZ"
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
  significance = '=IF(F2<0.0001,"****",IF(F2<0.001,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
  
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

##### ACTIVATION DURING RIPPLES ##### 
## load data file
fileN = "TableData_ripples_activation_placecells_RZ.txt"
dt = "rippleactivation_placecells_pairs_RZ"
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
  significance = '=IF(F2<0.0001,"****",IF(F2<0.001,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
  
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
