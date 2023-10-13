### load the libraries
# ALP 10/10/23
# This is to run the LME with interaction between frequencies (theta, slow gamma, medium gamma) and flicker condition
rm()

library(lme4)
library(emmeans)
library(sjPlot)
library(report)
library(ggplot2)

#load MATLAB output (in table format in .txt) file in R
foldN = "//ad.gatech.edu/bme/labs/singer/Abby/code/chronicflicker-ephys-prospectivecoding/results/LMM_R/"

##### PPC DATA! ##### 
## load data file
PPCtype <- c("VR")
celltype <- c("Pyr", "Int")
dtt <- c("ppc_crossreg_pyr", "ppc_crossreg_in")

for (PPCt in PPCtype) {
  for (CT in celltype) {
    fileN <- paste("TableData_PPC_", PPCt, "_CrossReg_", CT, ".txt", sep = "")
    dt <- paste("ppc_crossreg_", PPCt, "_", CT, sep = "")

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
    descript = 'PPC comparisons interaction frequency and group'
    comparison = "Flicker condition: 40 Hz vs. Random"
    testunit = "trials"
    significance = '=IF(F2<0.0001,"****",IF(F2<0.001,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel

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


