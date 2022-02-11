####################################################################################
##   Automatic Cleaning Pipeline for Luminex Data                                 ##
##                                                                                ##
##           Author: Dennis C.Y. Chan                                             ##
##           proctor.tools@gmail.com & dcc5251@psu.edu                            ##
##           GitHub: elizabethproctor                                             ##
##                                                                                ##
##   Please cite as:                                                              ##
##                                                                                ##
##                                                                                ##
##   This software is protected under the Gnu Public License, version 3 (GPLv3)   ##
####################################################################################

##Script for running logistic regression using nplr package. Purpose of this is to run logistic regression on reference
##multiplex data such that we have a regression relationship illustrating input fluorescence intensity to respective
##analyte concentration. 

gc() #Clean memory
rm(list=ls()) #Clean var space
Sys.setenv(JAVA_HOME='SET UP LOCATION OF JAVA MANUALLY HERE') #manually set Java location

#Required Libraries ----
library(nplr)
library(ggplot2)
library(gdata)
library(openxlsx)
library(sos)
library(xlsx)
library(MASS)
library(cluster)
library(glmnet)
library(factoextra)
library(goeveg)
library(dendextend)
library(nFactors)
library(philentropy)
#Additional Libraries ----
readLuminexExcelFileDir = '~/SET UP DIRECTORY TO VERSION OF readLimunexExcelFile_verX_0X.R'
source(readLuminexExcelFileDir)
#--Finish Set up of libraries--

#----Set up current directory and load file, and clean it----
wd = "SET UP WORKING DIRECTORY HERE" #Input working directory for target file to be loaded.
targetFileName = "file2" #Target File Name to load (no need for file type in the string variable)
beadTemp = readXLandSheet(wd, targetFileName, "Count")
netMFITemp = readXLandSheet(wd, targetFileName, 'Net MFI')
resultTemp = readXLandSheet(wd, targetFileName, 'Result')
classesTemp = readXLandSheet(wd, targetFileName, "Classes")
medianTemp = readXLandSheet(wd, targetFileName, "Median")

beadDataSet = beadTemp[[1]]
netMFIDataSet = netMFITemp[[1]]
resultDataSet = resultTemp[[1]]
classesDataSet = classesTemp[[1]]
medianDataSet = medianTemp[[1]]
refDataSet = medianDataSet 

#---- Removing temp values ----
remove(beadTemp)
remove(netMFITemp)
remove(resultTemp)
remove(classesTemp)
remove(medianTemp)

modeOfAggregation = "Mean" #Variable to select if Median/Mean is used to find the representative replicate value (options are as listed, includes first letter capital)
logFile = "textFileMean.log"
log_con = file(logFile) #Initializes and creates text log file

checkStandards = 0 #USER INPUT HERE TO CHECK/RECREATE STANDARDS FIRST

#---- Assuming that standards do not need to be checked since they were accepted by user from machine ----
if (checkStandards == 0){
  coefficientDataSet = readXLandSheet(wd, targetFileName, "Coefficients")[[1]]
  
  cat(paste("Date and Time of run: ", Sys.time(), '\n', sep=''), file=log_con)
  cat(paste("Running luminex file package Ver", sub('_', '.', strsplit(strsplit(readLuminexExcelFileDir, 'ver')[[1]][3], '[.]')[[1]][1]), "\n", sep=''), file=logFile, append=TRUE)
  cat("Assumed standards from multiplex output csv file is valid and results are also valid as a result from using the standard curves... \n", file=logFile, append=TRUE)
  cat("Cleaning up result data set first... \n", file=logFile, append=TRUE)
  resultDataSet = cleanNonExperimentVals(refDataSet, sampleIndex = 2, logFile) #Removes reference values rows/cells and compresses dataset 
  cat("Cleaning up bead count data set first... \n", file=logFile, append=TRUE)
  beadDataSet = cleanNonExperimentVals(beadDataSet, sampleIndex = 2, logFile) #Removes reference values rows/cells and compresses dataset 
  
  cleanRefLabels = c("Age") #User input to determine which class labels to examine before removing extrememties based off that label
  temp = cleanExtremeties(resultDataSet, classesDataSet, cleanRefLabels, 60, logFile) #Retrieves index where extremeties to be removed are stored
  if (any(is.na(temp)) == FALSE){
    resultDataSet[,temp] = NA #Sets those extremety columns with NA values for later preprocessing
  }
  resultDataSet = cleanNAN(resultDataSet) #Sets any N/A values and other iteratives of that to NA (type Numeric)
  beadDataSet = cleanNAN(beadDataSet) #Sets any N/A values and other iteratives of that to NA (type Numeric)
  
  cleanbeadDataSet = cleanCounts(beadDataSet, 20) #Examines bead count subset to set cells lower than threshold as NA
  mergedMatrixList = mergeSets(dataSet1 = resultDataSet,
                               dataSet2 = cleanbeadDataSet, 
                               naCol = temp, 
                               logFileName = logFile) #Merges cleaned counts with cleaned results to remove bad counts 
  resultDataSet = mergedMatrixList[[1]]
  merged_untaggedDataSet = mergedMatrixList[[2]]
  remove(mergedMatrixList)
  
  resultDataSet = resultDataSet[,-c(which(names(resultDataSet) == "Total.Events"))] #Removes total events column from matrix
  merged_untaggedDataSet = merged_untaggedDataSet[,-c(which(names(merged_untaggedDataSet) == "Total.Events"))]
  beadDataSet = beadDataSet[, names(resultDataSet)] #Sets beadDataSet to have the same columns as result Data Set. 
  
  resultDataSet = cbind(classesDataSet, resultDataSet, row.names = NULL)
  merged_untaggedDataSet = cbind(classesDataSet, merged_untaggedDataSet, row.names = NULL)
  beadDataSet = cbind(classesDataSet, beadDataSet, row.names = NULL)
  removeColVals = array(1:I(ncol(classesDataSet)+2)) #Appends the Location and Sample column from the result sheet to the classes sheet
  
  processedSet = cleanDataOffCVMat(dataSet = resultDataSet, 
                                       sampleIndex = length(removeColVals), 
                                       beadCount = beadDataSet,
                                       removeCols = removeColVals,
                                       tolerance = 25,  
                                       beadCountThresh = 20, 
                                       primerString = "CountLow",
                                       logFileName = logFile)
  
  #Creates the processed matrix which has the average/median calculated from the results. It is the totally clean data
  processedMat = calcAverageMat(processedDataSet = processedSet, 
                                sampleIndex = which(names(processedSet) == "Sample"), 
                                removeCols = removeColVals,
                                primer = c("CountLow", "All Low Count, Calc Avg", "Outlier"), 
                                mode=modeOfAggregation, 
                                logFileName = logFile)
  
  #Creates the untagged processed matrix which contains the averages calculated for observations and columns that had all low bead counts. 
  #Example would be say there exists a column for subject 001 where all replicates had low bead counts, an average/median is still calculated here,
  #this allows for retaining the subject/variable for analysis if desired. Thus, this helps create the 'dirty' data for exploratory stuff.
  merged_untaggedDataSet = formMatFromRefMatrix(dataSet = merged_untaggedDataSet, 
                                                refSet = processedSet, 
                                                primer = "All Low Count, Calc Avg")
  untaggedProcessedMat = calcAverageMat(processedDataSet = merged_untaggedDataSet, 
                                sampleIndex = which(names(merged_untaggedDataSet) == "Sample"), 
                                removeCols = removeColVals,
                                primer = c("CountLow", "All Low Count, Calc Avg", "Outlier"), 
                                mode=modeOfAggregation, 
                                logFileName = NA)
  
  #---- Writing the xlsx files for user to examine ----
  fileName = paste("Processed_Data_", targetFileName, "_", modeOfAggregation, ".xlsx", sep='')
  workbook = openxlsx::createWorkbook()
  addWorksheet(workbook, sheetName = 'ProcessedData')
  addWorksheet(workbook, sheetName = 'CleanedData_ObsOnly')
  addWorksheet(workbook, sheetName = 'CleanedData_VarsOnly')
  addWorksheet(workbook, sheetName = 'UncleanedData')
  writeData(workbook, sheet='ProcessedData',
            x=processedMat)
  writeData(workbook, sheet='CleanedData_ObsOnly', 
            x=processedMat[-unique(which(processedMat == "Low Bead Count Across All Replicates", arr.ind = TRUE)[,"row"]),])
  writeData(workbook, sheet='CleanedData_VarsOnly', 
            x=processedMat[,-unique(which(processedMat == "Low Bead Count Across All Replicates", arr.ind = TRUE)[,"col"])])
  writeData(workbook, sheet='UncleanedData', 
            x=untaggedProcessedMat)
  
  #Number of styles is equal to length of primer list as supposedly one style per primer string. 
  yellow_style = createStyle(fgFill="#FFFF00")
  red_style = createStyle(fgFill="#FF0000")
  orange_style = createStyle(fgFill="#FFB900")
  listOfStyles = list(yellow_style, red_style, orange_style)
  primerList = c("Outlier", "All Low Count, Calc Avg", "CountLow")
  listOfRemoveCols = c(1:which(names(processedSet) == "Sample"))
  tempTagMat = tagUntaggedMat(untaggedMat = processedSet,
                 compressMetric = "Sample",
                 removeCols = listOfRemoveCols,
                 primer = primerList)
  for (count in 1:length(listOfStyles)){
    tempMat = data.frame(which(tempTagMat == primerList[count], arr.ind = TRUE))
    y = unique(tempMat$col)
    for (i in 1:length(y)){
      x = tempMat[tempMat$col == y[i],1] # Extract only the row values 
      # x+1 as the first row is header, y+ncol(classesDataSet)+1 as ncol(classesDataSet) gives the class labels columns, and add 1 for the sample
      # column, do not account for the location column
      addStyle(workbook, sheet="UncleanedData", style=listOfStyles[[count]], rows=x+1, cols=y[i]+ncol(classesDataSet)+1, gridExpand=TRUE) 
    }
  }
  
  openxlsx::saveWorkbook(workbook, fileName, overwrite=TRUE)
  
  cat(paste("Created a file containing 4 different sheets. Filename is ", fileName, "\n", sep=''), file=logFile, append=TRUE)
  cat("Following are the sheet contents and description for the contents: \n", file=logFile, append=TRUE)
  cat(" -> ProcessedData: Contains the processed data matrix. However, observations from samples where all replicates for a given
      variable that all unsatisfactory bead counts were labeled as such in string form. \n", file=logFile, append=TRUE)
  cat("\n -> CleanedData_ObsOnly: Sheet contains a matrix that had removed all observations that had any instance of unsatisfactory 
      bead counts across all replicates for any specific variable. Example would be if subject X had a variable where all replicates for that
      variable had low bead counts. Thus subject X is removed from this matrix in this sheet. \n", file=logFile, append=TRUE)
  cat("\n -> CleanedData_VarsOnly: The same methodology from CleanedData_ObsOnly applies here with one key difference. Instead of removing the observation,
      the variables that contains instances of low bead counts across a subset of replicates is removed instead. \n", file=logFile, append=TRUE)
  cat("\n -> UncleanedData: Contains a data matrix that includes as many observations and columns as possible. In order to account for replicates
      that had unsatisfactory bead counts, an average/median was calculated from those replicates and used as the representative value to include
      in this data matrix. Thus, this is an uncleaned version of the processed matrix as there may be values contained here that were derived 
      from replicates that never had satisfcatory bead counts. Color Coding for UncleanedData sheet are as follows ", file=logFile, append=TRUE)
  cat("\n 
      ----> Red = All replicates had unsatisfactory bead counts \n
      ----> Orange = 1/2 Unsatisfactory Bead Counts in replicate subset, used the remaining satisfactory one \n
      ----> Yellow = Outlier found in replicate subset, used the remaining and calculated representative value \n", file=logFile, append=TRUE)
  cat("\n Finished! Cleared variables apart from resultDataSet and processedMat variables \n", file=logFile, append=TRUE)
  
  remove(beadDataSet, classesDataSet, cleanbeadDataSet, coefficientDataSet, listOfStyles, merged_untaggedDataSet, netMFIDataSet,
         yellow_style, orange_style, red_style, tempMat, tempTagMat, untaggedProcessedMat, refDataSet, processedSet)
  close(log_con)
}