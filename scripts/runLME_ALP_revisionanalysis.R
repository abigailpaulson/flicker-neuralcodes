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
  
 
  ##### Time 2 Next RZ - D2R ##### 
  ## this data structure only includes correct trials and post data
  ## load data file
  fileN = "TableData_Time2NextRZ_perTrial_postDataOnly_CtrlSpeed.txt"
  dt = "time2nextRZ"
  propTypes <-c("time_2_next_RZ")
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
    testData$animal = as.factor(testData$animal)
    testData$group = as.factor(testData$group)
    testData$Ctrl_speed_subset= as.factor(testData$Ctrl_speed_subset)
    
    #helpful info for stats (n per group)
    nUnit = c(sum(testData$Ctrl_speed_subset == 1), sum(testData$Ctrl_speed_subset == 2)) 
    
    ## loop over property types 
    ## set up helpful info for the model
    descript = 'time 2 next RZ per trial'
    comparison = "Speed subset: low vs. high"
    testunit = "trials"
    significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
    
    ## model, anova
    obj.lmer1 = lmerTest::lmer(time_2_next_RZ ~ Ctrl_speed_subset + (1|animal), data = testData)
    lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
    
    ## set up strings for the output
    `Group Size` = paste(gt, " -- speed 1 n =", as.character(nUnit[1]), "; speed 2 =", as.character(nUnit[2]))
    analysisUnit = "trials"
    `P-value` = round(lmer1anova$`Pr(>F)`, digits = 20)
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
  
  ##### THETA SEQUENCES - 180 deg PCR by Median Lick Number ##### 
  ## this data structure only includes correct trials and post data
  ## load data file
  fileN = "TableData_RRZPCR_180_perTrial_postDataOnly_LickNumber.txt"
  dt = "RRZPCR_median_LickNum"
  propTypes <-c("PCR_loc_trial")
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
    testData$animal = as.factor(testData$animal)
    testData$group = as.factor(testData$group)
    testData$Lick_number_subset= as.factor(testData$Lick_number_subset)
    
    #helpful info for stats (n per group)
    nUnit = c(sum(testData$Lick_number_subset == 1), sum(testData$Lick_number_subset == 2)) 
    
    ## loop over property types 
    ## set up helpful info for the model
    descript = 'time 2 next RZ per trial'
    comparison = "Speed subset: low vs. high"
    testunit = "trials"
    significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
    
    ## model, anova
    obj.lmer1 = lmerTest::lmer(PCR_loc_trial ~ Lick_number_subset + (1|animal), data = testData)
    lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
    
    ## set up strings for the output
    `Group Size` = paste(gt, " -- speed 1 n =", as.character(nUnit[1]), "; speed 2 =", as.character(nUnit[2]))
    analysisUnit = "trials"
    `P-value` = round(lmer1anova$`Pr(>F)`, digits = 20)
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

  
  ##### THETA SEQUENCES - Per Trial - Combined Groups - Unengaged vs. Correct Trials ##### 
  
  ## load data file
  fileN = "TableData_TrialData_180Decoding_ThetaSeq_Unengaged_Correct_Trials.txt"
  dt = "thetaseq_Unengaged_correct_pertrial"
  propTypes <-c("PCR_loc_trial","Ctrl_speed")
  
  #initialize table structure to add stats to
  outputStats = data.frame()
  
  tblFilename = paste0(foldN, fileN)
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create text file to save LMM outputs 
  csvFilename = paste0(foldN, "/stats_lme4_", dt, '.csv') #create a csv that can be used for supplement table
  
  #read data
  mydata <- read.csv(tblFilename, head=T)
  
  ## get data subset of interest - this struct includes only plotted data, so include everything
  testData <- mydata
  
  #factorize variables to be treated as categories/levels instead of numerical value
  testData$animal = as.factor(testData$animal)
  testData$group = as.factor(testData$group)
  testData$isUnengaged = as.factor(testData$isUnengaged) # if not unengaged then it is correct
  
  #helpful info for stats (n per group)
  nUnit = c(sum(testData$isUnengaged == 1), sum(testData$isUnengaged == 0)) 
  
  ## loop over property types 
  ## set up helpful info for the model
  for (pt in propTypes) {
    descript = pt
    comparison = "unengaged vs. correct trials both groups "
    testunit = "trial"
    significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
    
    ## model, anova
    obj.lmer1 = lmerTest::lmer(get(pt) ~ isUnengaged + (1|animal), data = testData)
    lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
    
    ## set up strings for the output
    `Group Size` = paste("unengaged n =", as.character(nUnit[1]), "; correct n =", as.character(nUnit[2]))
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
  
  
  ##### THETA SEQUENCES - Per Trial - By Group - Unengaged vs. Correct Trials ##### 
  ## this data structure only includes correct trials and post data
  ## load data file
  fileN = "TableData_TrialData_180Decoding_ThetaSeq_Unengaged_Correct_Trials.txt"
  dt = "thetaseq_Unengaged_correct_pertrial_pergroup"
  propTypes <-c("PCR_loc_trial")
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
    
    #factorize variables to be treated as categories/levels instead of numerical value
    testData$animal = as.factor(testData$animal)
    testData$group = as.factor(testData$group)
    testData$isUnengaged = as.factor(testData$isUnengaged)
    
    #helpful info for stats (n per group)
    nUnit = c(sum(testData$isUnengaged == 1), sum(testData$isUnengaged == 0)) 
    
    ## loop over property types 
    ## set up helpful info for the model
    descript = 'time 2 next RZ per trial'
    comparison = "Speed subset: low vs. high"
    testunit = "trials"
    significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
    
    ## model, anova
    obj.lmer1 = lmerTest::lmer(PCR_loc_trial ~ isUnengaged + (1|animal), data = testData)
    lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
    
    ## set up strings for the output
    `Group Size` = paste(gt, " -- unengaged n =", as.character(nUnit[1]), "; engaged =", as.character(nUnit[2]))
    analysisUnit = "trials"
    `P-value` = round(lmer1anova$`Pr(>F)`, digits = 20)
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
  
  
  
  ##### THETA SEQUENCES - Per Day - By Group - Unengaged vs. Correct Trial Difference ##### 
  ## this data structure only includes correct trials and post data
  ## load data file
  rm(mydata)
  rm(testData)
  rm(obj.lmer1)
  rm(lmer1anova)
  fileN = "TableData_DayData_180Decoding_ThetaSeq_Unengaged_Correct.txt"
  dt = "thetaseq_Unengaged_correct_perDay_pergroup"
  propTypes <-c("avg_PCR_RRZ_difference")
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
    testData <- mydata
    
    #factorize variables to be treated as categories/levels instead of numerical value
    testData$animal = as.factor(testData$animal)
    testData$group = as.factor(testData$group)

    
    #helpful info for stats (n per group)
    nUnit = c(sum(testData$group == "random"), sum(testData$group == "gamma")) 
    
    ## loop over property types 
    ## set up helpful info for the model
    descript = 'time 2 next RZ per trial'
    comparison = "Speed subset: low vs. high"
    testunit = "trials"
    significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
    
    ## model, anova
    obj.lmer1 = lmerTest::lmer(avg_PCR_RRZ_difference ~ group + (1|animal), data = testData)
    lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
    
    ## set up strings for the output
    `Group Size` = paste(" -- random n =", as.character(nUnit[1]), "; 40 Hz =", as.character(nUnit[2]))
    analysisUnit = "trials"
    `P-value` = round(lmer1anova$`Pr(>F)`, digits = 20)
    `F-value` = paste("F", as.character(round(lmer1anova$NumDF, digits=5)), ",", as.character(round(lmer1anova$DenDF, digits = 5)), " = ", as.character(round(lmer1anova$`F value`, digits = 5)), sep="")
    
    ## compile into output frame
    tmpStats = data.frame(descript, testunit, comparison, `Group Size`, `P-value`, significance,`F-value`)
    outputStats = rbind(outputStats, tmpStats)
    
    ## output all stats details
    statStr = paste(gt, "|", descript, "|", comparison, "|", "n =", as.character(nUnit[1]), ", n =", as.character(nUnit[2]), testunit, "|", "p =", as.character(round(lmer1anova$`Pr(>F)`, digits = 5)), "|", "F", as.character(round(lmer1anova$NumDF, digits = 5)), as.character(round(lmer1anova$DenDF, digits = 5)), "=", as.character(round(lmer1anova$`F value`, digits = 5)))
    
    sink(txtFilename, append = TRUE)
    print(report(obj.lmer1))
    print(statStr)
    print(summary(obj.lmer1))
    print(anova(obj.lmer1, ddf="Kenward-Roger"))
    print ("-------------------------------------------------------")
    sink()
    closeAllConnections()
    
    write.csv(outputStats, csvFilename)
    

  
  
    
    ##### THETA SEQUENCES - Per Trial - By Group - RRZ vs Ctrl Region ##### 
    ## this data structure only includes correct trials and post data
    ## load data file
    fileN = "TableData_TrialData_180Decoding_ThetaSeq_RRZ_vs_Control.txt"
    dt = "thetaseq_RRZ_vs_Ctrl_pertrial_pergroup"
    propTypes <-c("PCR_trial")
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
      
      #factorize variables to be treated as categories/levels instead of numerical value
      testData$animal = as.factor(testData$animal)
      testData$group = as.factor(testData$group)
      testData$PCRtype = as.factor(testData$PCRtype)
      
      #helpful info for stats (n per group)
      nUnit = c(sum(testData$group == gt)/2, sum(testData$group == gt)/2) 
      
      ## loop over property types 
      ## set up helpful info for the model
      descript = 'PCR RRZ vs Ctrl'
      comparison = "RRZ vs. Ctrl"
      testunit = "trials"
      significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
      
      ## model, anova
      obj.lmer1 = lmerTest::lmer(PCR_trial ~ PCRtype + (1|animal), data = testData)
      lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
      
      ## set up strings for the output
      `Group Size` = paste(gt, " -- n =", as.character(nUnit[1]), "; n =", as.character(nUnit[2]))
      analysisUnit = "trials"
      `P-value` = round(lmer1anova$`Pr(>F)`, digits = 20)
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
    
    ##### THETA SEQUENCES - Per Trial - By Group - Unoccupied zone decoding vs. not - Ctrl Speed ##### 
    ## this data structure only includes correct trials and post data
    ## load data file
    fileN = "TableData_TrialData_360Decoding_ThetaSeq_PCR_OtherZoneDecoding.txt"
    dt = "thetaseq_RRZ_OtherZoneDecoding_pertrial_pergroup_CtrlSpeed"
    propTypes <-c("Ctrl_speed")
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
      
      #factorize variables to be treated as categories/levels instead of numerical value
      testData$animal = as.factor(testData$animal)
      testData$group = as.factor(testData$group)
      testData$sig_otherpos_RRZ= as.factor(testData$sig_otherpos_RRZ)
      
      #helpful info for stats (n per group)
      nUnit = c(sum(testData$group == gt), sum(testData$group == gt)) 
      
      ## loop over property types 
      ## set up helpful info for the model
      descript = 'Ctrl speed Other zone vs. not'
      comparison = "Other zone decoding vs not"
      testunit = "trials"
      significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
      
      ## model, anova
      obj.lmer1 = lmerTest::lmer(Ctrl_speed ~ sig_otherpos_RRZ + (1|animal), data = testData)
      lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
      
      ## set up strings for the output
      `Group Size` = paste(gt, " -- n =", as.character(nUnit[1]), "; n =", as.character(nUnit[2]))
      analysisUnit = "trials"
      `P-value` = round(lmer1anova$`Pr(>F)`, digits = 20)
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

    ##### THETA SEQUENCES - Per Trial - By Group - Unoccupied zone decoding vs. not -lickDI ##### 
    ## this data structure only includes correct trials and post data
    ## load data file
    fileN = "TableData_TrialData_360Decoding_ThetaSeq_PCR_OtherZoneDecoding.txt"
    dt = "thetaseq_RRZ_OtherZoneDecoding_pertrial_pergroup_LickDI"
    propTypes <-c("LickDI")
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
      
      #factorize variables to be treated as categories/levels instead of numerical value
      testData$animal = as.factor(testData$animal)
      testData$group = as.factor(testData$group)
      testData$sig_otherpos_RRZ= as.factor(testData$sig_otherpos_RRZ)
      
      #helpful info for stats (n per group)
      nUnit = c(sum(testData$group == gt), sum(testData$group == gt)) 
      
      ## loop over property types 
      ## set up helpful info for the model
      descript = 'LIckDI Other zone vs. not'
      comparison = "Other zone decoding vs not"
      testunit = "trials"
      significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
      
      ## model, anova
      obj.lmer1 = lmerTest::lmer(LickDI ~ sig_otherpos_RRZ + (1|animal), data = testData)
      lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
      
      ## set up strings for the output
      `Group Size` = paste(gt, " -- n =", as.character(nUnit[1]), "; n =", as.character(nUnit[2]))
      analysisUnit = "trials"
      `P-value` = round(lmer1anova$`Pr(>F)`, digits = 20)
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
    
    ##### THETA SEQUENCES - Per Trial - By Group - Unoccupied zone decoding vs. not - lick latency ##### 
    ## this data structure only includes correct trials and post data
    ## load data file
    fileN = "TableData_TrialData_360Decoding_ThetaSeq_PCR_OtherZoneDecoding.txt"
    dt = "thetaseq_RRZ_OtherZoneDecoding_pertrial_pergroup_LickLatency"
    propTypes <-c("licklatency_s")
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
      
      #factorize variables to be treated as categories/levels instead of numerical value
      testData$animal = as.factor(testData$animal)
      testData$group = as.factor(testData$group)
      testData$sig_otherpos_RRZ= as.factor(testData$sig_otherpos_RRZ)
      
      #helpful info for stats (n per group)
      nUnit = c(sum(testData$group == gt), sum(testData$group == gt)) 
      
      ## loop over property types 
      ## set up helpful info for the model
      descript = 'lick latency Other zone vs. not'
      comparison = "Other zone decoding vs not"
      testunit = "trials"
      significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
      
      ## model, anova
      obj.lmer1 = lmerTest::lmer(licklatency_s ~ sig_otherpos_RRZ + (1|animal), data = testData)
      lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
      
      ## set up strings for the output
      `Group Size` = paste(gt, " -- n =", as.character(nUnit[1]), "; n =", as.character(nUnit[2]))
      analysisUnit = "trials"
      `P-value` = round(lmer1anova$`Pr(>F)`, digits = 20)
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
    
    ##### Ripples - Per Trial - By Group - Unoccupied zone decoding vs. not - Ctrl Speed ##### 
    ## this data structure only includes correct trials and post data
    ## load data file
    fileN = "TableData_TrialData_360Decoding_Ripples_otherzone_pertrial.txt"
    dt = "ripples_RZ_OtherZoneDecoding_pertrial_pergroup_CtrlSpeed"
    propTypes <-c("Ctrl_speed")
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
      
      #factorize variables to be treated as categories/levels instead of numerical value
      testData$animal = as.factor(testData$animal)
      testData$group = as.factor(testData$group)
      testData$isNonLocal = as.factor(testData$isNonLocal)
      
      #helpful info for stats (n per group)
      nUnit = c(sum(testData$isNonLocal == 0), sum(testData$isNonLocal == 1)) 
      
      ## loop over property types 
      ## set up helpful info for the model
      descript = 'Ctrl speed Other zone vs. not'
      comparison = "Other zone decoding vs not"
      testunit = "trials"
      significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
      
      ## model, anova
      obj.lmer1 = lmerTest::lmer(Ctrl_speed ~ isNonLocal + (1|animal), data = testData)
      lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
      
      ## set up strings for the output
      `Group Size` = paste(gt, " -- n =", as.character(nUnit[1]), "; n =", as.character(nUnit[2]))
      analysisUnit = "trials"
      `P-value` = round(lmer1anova$`Pr(>F)`, digits = 20)
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
    
    ##### Ripples - Per Trial - By Group - Unoccupied zone decoding vs. not - Lick DI ##### 
    ## this data structure only includes correct trials and post data
    ## load data file
    fileN = "TableData_TrialData_360Decoding_Ripples_otherzone_pertrial.txt"
    dt = "ripples_RZ_OtherZoneDecoding_pertrial_pergroup_lickDI"
    propTypes <-c("lickDI")
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
      
      #factorize variables to be treated as categories/levels instead of numerical value
      testData$animal = as.factor(testData$animal)
      testData$group = as.factor(testData$group)
      testData$isNonLocal = as.factor(testData$isNonLocal)
      
      #helpful info for stats (n per group)
      nUnit = c(sum(testData$isNonLocal == 0), sum(testData$isNonLocal == 1)) 
      
      ## loop over property types 
      ## set up helpful info for the model
      descript = 'Lick DI Other zone vs. not'
      comparison = "Other zone decoding vs not"
      testunit = "trials"
      significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
      
      ## model, anova
      obj.lmer1 = lmerTest::lmer(lickDI ~ isNonLocal + (1|animal), data = testData)
      lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
      
      ## set up strings for the output
      `Group Size` = paste(gt, " -- n =", as.character(nUnit[1]), "; n =", as.character(nUnit[2]))
      analysisUnit = "trials"
      `P-value` = round(lmer1anova$`Pr(>F)`, digits = 20)
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
    
    ##### Ripples - Per Trial - By Group - Unoccupied zone decoding vs. not - licklatency ##### 
    ## this data structure only includes correct trials and post data
    ## load data file
    fileN = "TableData_TrialData_360Decoding_Ripples_otherzone_pertrial.txt"
    dt = "ripples_RZ_OtherZoneDecoding_pertrial_pergroup_licklatency"
    propTypes <-c("licklatency_s")
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
      
      #factorize variables to be treated as categories/levels instead of numerical value
      testData$animal = as.factor(testData$animal)
      testData$group = as.factor(testData$group)
      testData$isNonLocal = as.factor(testData$isNonLocal)
      
      #helpful info for stats (n per group)
      nUnit = c(sum(testData$isNonLocal == 0), sum(testData$isNonLocal == 1)) 
      
      ## loop over property types 
      ## set up helpful info for the model
      descript = 'Lick Latency Other zone vs. not'
      comparison = "Other zone decoding vs not"
      testunit = "trials"
      significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
      
      ## model, anova
      obj.lmer1 = lmerTest::lmer(licklatency_s ~ isNonLocal + (1|animal), data = testData)
      lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
      
      ## set up strings for the output
      `Group Size` = paste(gt, " -- n =", as.character(nUnit[1]), "; n =", as.character(nUnit[2]))
      analysisUnit = "trials"
      `P-value` = round(lmer1anova$`Pr(>F)`, digits = 20)
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
    
    ##### proportion of trials with unoccupied zone decoding by fast/slow ##### 
    ## this data structure only includes correct trials and post data
    ## load data file
    fileN = "TableData_DayData_360Decoding_ThetaSeq_PCR_OtherZoneDecoding.txt"
    dt = "thetaseq_RRZ_OtherZoneDecoding_perday_pergroup_splitbytrialspeed"
    propTypes <-c("prop_otherZone_byspeed")
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
      
      #factorize variables to be treated as categories/levels instead of numerical value
      testData$animal = as.factor(testData$animal)
      testData$group = as.factor(testData$group)
      testData$speedType = as.factor(testData$speedType)
      
      #helpful info for stats (n per group)
      nUnit = c(sum(testData$speedType == 0), sum(testData$speedType == 1)) 
      
      ## loop over property types 
      ## set up helpful info for the model
      descript = 'prop trials with unoccupied zone decoding by speed type'
      comparison = "within group across speed types"
      testunit = "days"
      significance = '=IF(F2<0.001,"****",IF(F2<0.005,"***",IF(F2<0.01,"**",IF(F2<0.05,"*","ns"))))' #this is the formula for the excel file, it will evaluate in excel
      
      ## model, anova
      obj.lmer1 = lmerTest::lmer(propNonLocal ~ speedType + (1|animal), data = testData)
      lmer1anova = anova(obj.lmer1, ddf="Kenward-Roger")
      
      ## set up strings for the output
      `Group Size` = paste(gt, " -- n =", as.character(nUnit[1]), "; n =", as.character(nUnit[2]))
      analysisUnit = "trials"
      `P-value` = round(lmer1anova$`Pr(>F)`, digits = 20)
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