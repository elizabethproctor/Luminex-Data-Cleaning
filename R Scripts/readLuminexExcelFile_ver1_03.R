####################################################################################
##   'Package' that holds in house R functions for automatic cleaning pipeline    ##
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
#R Script that contains the functions that read in a luminex excel file that has already been segmented into sheets containing specific data types
#It contains the scripts needed to clean the data if required, allows for finding columns with zero variance etc. 
#
#----Version 1.01 had the following functions:----
# + is.interger0
# + readXLandSheet
# + cellColor
# + cleanExcelFile
# + decimalPlaces
# + addTolerance
# + removeZeroVariance
# + nCr
# + findDiff
# + cleanDataOffCVMat
# + selectionVals
# + locateOutlier
# + findAllPairs
# + extractRelevantVals
# + createCVMat
# + calculateCoeffVar
# + cleanCounts
# + cleanExtremeties
# + cleanNAN
# + mergeSets
# + cleanNonExperimentVals
# + formMatFromRefMatrix
# + compressMatrix
# + getMode
# + calcAverageMat
# + tagUntaggedMat
#
#----Version 1.02 had the following changes made:----
# + added magnitude_vector
# + added z_score_standardization
# + added Libraries that were required.

#----Version 1.03 had the following changes made:----
# + added lines at 591-596 to trim white space at beginning of every cell within Sample column, as leading white space causes
# the numeric represnetation of the cell to be negative, and causes the entire column to be removed, which is conterproductive to
# the functions that depend on mergeSet. 

#Version 1.03 Wrapped at 2-3-2021

#---- Adding Libraries ----
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
library(matrixStats)
library(philentropy)
library(matlib)

#---- Functions ----
#Checks for interger0 type
is.interger0 = function(inputVar){
  is.integer(inputVar) && !length(inputVar)
}

#Reads xlsx file and determines which sheet to read dependent on inputs
readXLandSheet = function(wd, targetFileName, sheetTargetName){
  setwd(wd)
  currentFiles = list.files()
  for (i in 1:length(currentFiles)){
    localFile = tools::file_path_sans_ext(currentFiles[i])
    if (localFile == targetFileName){
      sheetName = getSheetNames(currentFiles[i])
      for (j in 1:length(sheetName)){
        if (sheetName[j] == sheetTargetName){
          dataSet = openxlsx::read.xlsx(currentFiles[i], j)
          workbookObj = openxlsx::loadWorkbook(currentFiles[i])
          sheet = openxlsx::readWorkbook(workbookObj, sheet = sheetTargetName)
          output = list(dataSet, sheet)
          return(output)
        }
      }
    } 
  }
  print("Working Directory does not have file, quitting script... Please reset directory")
  return(NA)
}

#Examines cell color from excel file
cellColor = function(style) {
  fg  = style$getFillForegroundXSSFColor()
  rgb = tryCatch(fg$getRgb(), error = function(e) NULL)
  rgb = paste(rgb, collapse = "")
  return(rgb)
}

#Cleans input excel file
cleanExcelFile = function(sheet, dataSet){
  rows = getRows(sheet[[1]])
  cells = getCells(rows)
  styles = sapply(cells, getCellStyle)
  cellStyleVals = sapply(styles, cellColor)
  pheno = list(normal = "", remove = "ff0000") #HIGHLIGHT RED IN EXCEL
  m = match(sapply(styles, cellColor), pheno)
  dataNames = names(pheno)[m]
  tempMatrix = matrix(NA, nrow = nrow(dataSet), ncol=ncol(dataSet))
  tempX = 1
  tempY = 1
  for (i in I(ncol(dataSet)+1):length(m)){
    continue = F
    while (continue == F){
      if (is.na(dataSet[tempX, tempY]) == F && dataNames[i] == "remove"){
        dataSet[tempX, tempY] = NA
        continue = T
      } else if (is.na(dataSet[tempX, tempY]) == T){
        tempY = tempY + 1
      } else if (is.na(dataSet[tempX, tempY]) == F && dataNames[i] == "normal"){
        continue = T
      }
      if (tempY > ncol(dataSet)){
        tempY = 1
        tempX = tempX + 1
      }
      if (tempX > nrow(dataSet)){break}
    }
    tempY = tempY + 1
    if (tempY > ncol(dataSet)){
      tempY = 1
      tempX = tempX + 1
    }
  }
  return(dataSet)
}

#Checks number of decimal places
decimalPlaces = function(number){
  if (abs(number - round(number)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(number)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

#Creates a number greater by one unit above the number of decimal places
addTolerance = function(number){
  output = '0.'
  for (i in 1:I(number-1)){
    output = paste(output, '0', sep='')
  }
  output = paste(output, '1', sep='')
  return(as.numeric(as.character(output)))
}

#Outputs zero variance indices to remove 
removeZeroVariance = function(initialRemoveIndex, subset){
  removeIndex = matrix(NA, ncol(subset)) #Initialization of array space
  count = 1
  for (i in 1:length(initialRemoveIndex)){
    removeIndex[count] = initialRemoveIndex[i]
    count = count + 1
  }
  colI_Index = c(1:ncol(subset))
  colI_Index = colI_Index[-(initialRemoveIndex)]
  for (colI in colI_Index){
    if (var(subset[!is.na(subset[[colI]]),colI]) == 0){
      removeIndex[count] = as.character(colI)
      count = count + 1
    }
  }
  removeIndex = removeIndex[!is.na(removeIndex)]
  removeIndex = as.numeric(removeIndex)
  return(removeIndex)
}

#Does combination formula (nCr)
nCr = function(n, r){
  output = factorial(n)/I(factorial(r) * factorial(n-r))
  return(output)
}

#Finds the maximum/minimum difference of elements within array
findDiff = function(array, option, removeNA){
  if (removeNA == T){
    array = array[!is.na(array)]
  }
  array = as.numeric(array[-grep("CountLow", array)])
  meanRef = mean(as.numeric(array))
  tempArray = array(NaN, length(array))
  for (i in 1:length(array)){
    tempArray[i] = abs(meanRef - array[i])
  }
  if(option == "max"){
    output = which.max(tempArray)
  } else if (option == "min"){
    output = which.min(tempArray)
  }
  return(output)
}

#Inputs CV matrix and compares with full rank matrix to clean data by removing ones
#that deviate the most, objective is to minimize the cV based off concentrations
cleanDataOffCVMat = function(dataSet, sampleIndex, beadCount, removeCols, tolerance, beadCountThresh, primerString, logFileName){
  #cat("\n#---- Calcualting Coefficient of variance matrices on input matrix ----# \n", file = logFileName, append = T)
  
  cleanedSet = dataSet
  cvMat = createCVMat(cleanedSet, sampleIndex, removeCols, primerString)
  
  # intraAssayCV = array(NA, ncol(cvMat))
  # for (i in 2:length(intraAssayCV)){
  #   intraAssayCV[i] = mean(cvMat[,i], na.rm=T)
  # }
  # names(intraAssayCV) = names(cvMat)
  # intraAssayCV = intraAssayCV[-is.na(intraAssayCV)]
  
  for (i in 1:nrow(cvMat)){
    currentSubject = cvMat[i,][[1]] #First column of CV will always be sample ID. 
    currentMat = cleanedSet[cleanedSet[[sampleIndex]] == currentSubject, ] #Creates a temp matrix for the current subject/class only.
    currentBead = beadCount[beadCount[[sampleIndex]] == currentSubject, ]
    currentCVMat = cvMat[i,]
    if (any(is.na(currentCVMat)) == F){
      countLowArrayIndex = array(NA, nrow(currentMat))
      for (j in 1:length(countLowArrayIndex)){
        numberOfLowCounts = length(grep(primerString, currentMat[j,]))
        countLowArrayIndex[j] = numberOfLowCounts
      }
      currentMat = selectionVals(currentCVMat, currentMat, currentBead, countLowArrayIndex, removeCols, beadCountThresh, threshold=tolerance)
    } else { #There is a technical replicate set that have low bead accounts across all observations 
      currentMat = selectionVals(currentCVMat, currentMat, currentBead, countLowArrayIndex, removeCols, beadCountThresh, threshold=tolerance)
    }
    
    cleanedSet[cleanedSet$Sample == currentSubject, -removeCols] = currentMat
  }
  return(cleanedSet)
}

#Function that takes in the cv matrix, current matrix, current bead matrix and an array indicating how many countLow designations per
#observation. 
selectionVals = function(currentCVMat, currentMat, currentBead, countLowArrayIndex, removeCols, beadCountThresh, threshold){
  minIndex =  which(countLowArrayIndex == countLowArrayIndex[which.min(countLowArrayIndex)])
  labelCols = currentMat[, removeCols]
  currentMat = currentMat[, -removeCols]
  currentBead = currentBead[, -removeCols]
  currentCVMat = currentCVMat[, -1]  #First column is sample value
  for (i in 1:ncol(currentCVMat)){
    if (is.na(currentCVMat[i]) == F){
      if(currentCVMat[i] <= threshold){
        #Do nothing since the current CV for the variable is good, simply log it (NEED TO WRITE CODE FOR THIS) ----
        # if (is.interger0(grep("CountLow", currentMat[,i])) == F){
        #   currentMat[grep("CountLow", currentMat[,i]), i] = "Remove"
        # } else { } 
        currentMat[,i] = currentMat[,i] #Do nothing ----
      } else {
        badCountRow = which("CountLow" == currentMat[,i]) #Retrieves row index if it has countLOw index
        if (is.interger0(badCountRow) == T){ #So no cells with bad bead counts
          currentMat[,i] = locateOutlier(currentMat[,i], 100)
        } else if (length(badCountRow) >= 3){ #If all rows have bad bead counts then the following statement will be executed
          currentMat[,i] = "Remove"
        } else { #Everything else is left alone since the next function is calcAverage which does averages with existing values not as strings
          currentMat[,i] = currentMat[,i]
        }
      }
    } else {
      #Append log file detailing which one is NA index
      NAIndex = which(is.na(currentCVMat))
      currentMat[,i] = "All Low Count, Calc Avg"
    }
  }
  return(currentMat)
  
  #---- Draft ----
  # if (any(countLowArrayIndex > 0)){
  #   minIndex =  which(countLowArrayIndex == countLowArrayIndex[which.min(countLowArrayIndex)])
  #   if (length(minIndex) > 2){
  #     currentMat = extractRelevantVals(currentCVMat, currentMat, currentBead, primer = "CountLow", removeCols, beadCountThresh, threshold)
  #   } else if (length(minIndex) == 1){
  #     return(currentMat[minIndex,]) #Return the observation with the lowest number of bad bead counts
  #   } else if (length(minIndex) == 2){
  #     currentMat = extractRelevantVals(currentCVMat, currentMat[minIndex,], currentBead[minIndex,], primer = "CountLow", removeCols, beadCountThresh, threshold)
  #   } else if (length(minIndex) == 0){
  #     print("ERROR") #No minimum value is an error and bug. 
  #   }
  # } else {
  #   return(currentMat) #Do nothing and return base currentMat
  # }
}

#Function that allows to investigate a given vector and find pair-wise differences (absolute differences). The method to find an outlier
#would be to use these pair-wise differences and find if they differ significantly in magnitude, and if so, how many pair-wise differences
#fall into that category, and if there is a common denominator in this subset. The common-denominator is then the outlier by definition.
locateOutlier = function(colMat, threshold){
  numOfPairs = nCr(length(colMat), 2)
  tempArray = array(NA, numOfPairs)
  names(tempArray) = findAllPairs(1, length(colMat), interval=1)
  for (i in 1:numOfPairs){
    valA = colMat[as.numeric(as.character((strsplit(names(tempArray), ',')[[i]][1])))]
    valB = colMat[as.numeric(as.character((strsplit(names(tempArray), ',')[[i]][2])))]
    tempArray[i] = abs(as.numeric(valA) - as.numeric(valB))
    if (tempArray[i] < 1){tempArray[i] = 1} #Done for easier mathematical calcuations and debugging
  }
  
  #There are three conditions that can occur with 3 values that arose from pair-wise subtraction. 
    #Either they are all equally large/small
    #One small and two large
    #Or all random and hold no pattern
  #Therefore, the only discernable pattern to look for is where there is one significantly lower than the other two.
  lowestDiff = which(tempArray == min(tempArray))
  minVal = tempArray[lowestDiff]
  tempArray = tempArray[-lowestDiff]
  ratioA = as.numeric(unname(tempArray[1]))/as.numeric(unname(minVal))
  ratioB = as.numeric(unname(tempArray[2]))/as.numeric(unname(minVal))
  
  if(is.na(ratioA) == TRUE || is.na(ratioB) == TRUE){
    print("y")
  }
  
  if ( (any( (as.numeric(colMat) < 100) == FALSE)) == FALSE){
    tempOutput = colMat #Do nothing
  } else if (ratioA > 2 & ratioB > 2){ #If the ratios of the two differences over min difference is greater than 2
    excludeIndex = getmode(as.numeric(strsplit(paste(names(tempArray[1]), names(tempArray[2]), sep=','), ',')[[1]]))
    tempOutput = colMat 
    tempOutput[excludeIndex] = "Outlier"
  } else {
    #Otherwise do nothing since the differences do not follow a specific pattern and it is possible the points are simply spread out. 
    tempOutput = colMat
  }
  return(tempOutput)
}

#Subfunction that is most likely utilized in locateOutlier, it finds all possible pairs given an input starting index, a final index, and 
#an interval variable that helps determine if there is any numbers that are to be skipped etc. Returns pairs in string format
findAllPairs = function(startNum, endNum, interval){
  allNums = seq(startNum, endNum, interval)
  numOfPairs= as.numeric(nCr(length(allNums), 2))
  tempOutput = array(NA, numOfPairs)
  i = 1
  for (j in 1:length(allNums)){
    firstNumOfPair = as.character(allNums[j])
    for (k in I(j+1):length(allNums)){
      secondNumOfPair = as.character(allNums[k])
      tempString = paste(firstNumOfPair, secondNumOfPair, sep=',')
      tempOutput[i] = tempString
      i = i + 1
    }
    if(any(is.na(tempOutput)) == F){break}
  }
  return(tempOutput)
}

#Function that intakes the coefficeint of variance matrix, the current matrix of concentrations and the current bead count for that sample
#Objective is to narrow down the curMat input by examining the bead counts and applying the values to the following conditions:
#- First inspect if the CV crosses the threshold value, if so determine how many good bead counts there are for that variable 
#- Find if the low bead counts exist in the same variable
#- If so, return all observations 
#- If not, determine what the differences are for the bead counts for each observation from the threshold and score them, do this only
#  with variables that have low bead counts but are not present across all observations. 
extractRelevantVals = function(cvMat, curMat, curBead, primer, removeCols, beadCountThresh, threshold){
  classLabels = curMat[,removeCols]
  curMat = curMat[,-removeCols]
  curBead = curBead[,-removeCols]
  countLowInstances_inColumn = colSums(curMat == primer) #Counts number of CountLow tags per variable
  maxVal = max(countLowInstances_inColumn) #Finds maximum number of CountLow tags per variable
  minVal = min(countLowInstances_inColumn[countLowInstances_inColumn > 0]) #Finds minimum number of CountLow tags per variable
  if (minVal == maxVal){ #Then the low bead counts are in the same variables, simple case, return them all
    return(curMat)
  } else { #Low Bead counts are not in the same variables, requires further processing
    tempArray = array(NA, length(countLowInstances_inColumn[countLowInstances_inColumn > 0])) #Sets array length for number of variables
    names(tempArray) = names(countLowInstances_inColumn[countLowInstances_inColumn > 0])
    for (j in 1:length(tempArray)){
      tempArray[j] = which(curBead[, names(tempArray)[j]] == max(curBead[, names(tempArray)[j]]))
    }
    desiredRow = getmode(tempArray)
    curMat = cbind(classLabels, curMat)
    curMat = curMat[desiredRow,]
    return(curMat)
  }
}

#Creates coefficient of variation (CV) matrix 
createCVMat = function(processedDataSet, sampleIndex, removeCols, primerString){
  subjectVal = unique(processedDataSet[[sampleIndex]])
  cvRawMat = matrix(NA, nrow = length(subjectVal), ncol = ncol(processedDataSet))
  cvRawMat = as.data.frame(cvRawMat)
  colnames(cvRawMat) = colnames(processedDataSet)
  labelCols = unique(processedDataSet[, c(sampleIndex)])
  cvRawMat = cvRawMat[, -removeCols]
  for (a in 1:length(subjectVal)){
    currentSubject = subjectVal[a]
    tempMatrix = processedDataSet[processedDataSet[[sampleIndex]] == currentSubject,]
    tempMatrix = tempMatrix[,-removeCols]
      if(nrow(tempMatrix) > 0){
        for (b in 1:ncol(tempMatrix)){
          cvMat = calculateCoeffVar(tempMatrix[,b], primerString)
          cvRawMat[a, b] = cvMat[3]
        }
      } else {
        print("ERROR IN CREATING CV MAT")
      }
  }
  cvRawMat = cbind(labelCols, cvRawMat) #Created a CV matrix for each individual result
  return(cvRawMat)
}

#Calculates the coefficient of variation. Accounts for string variables that are to be removed as well. Output is a 4 cell array
#well the first cell is the average, second is median, third cell is the standard deviation, and the fourth cell is the CV value (in percentage).
calculateCoeffVar = function(dataSet, primerString){
  cvReturn = array(NA, 4)
  #Loops through string array to look for any instances of each primer string in current dataset to retrieve and examine
  for (i in 1:length(primerString)){
    removedIndices = grep(primerString[i], dataSet)
    if (is.interger0(removedIndices) == FALSE){ break }
  }
  if (is.interger0(removedIndices) == FALSE){
    if (length(dataSet[-removedIndices]) > 1){
      average = mean(as.numeric(dataSet[-removedIndices]))
      median = median(as.numeric(dataSet[-removedIndices]))
      standardDeviation = sd(as.numeric(dataSet[-removedIndices]))
    } else {
      if (length(dataSet[-removedIndices]) == 1){
        average = mean(as.numeric(dataSet[-removedIndices]))
        median = median(as.numeric(dataSet[-removedIndices]))
        standardDeviation = 0
      } else {
        average = NA
        median = NA
        standardDeviation = NA
      }
    }
  } else { #Criteria is only fulfilled if removed indices is empty, meaning that there is no  count empty criteria 
    average = mean(as.numeric(dataSet))
    median = median(as.numeric(dataSet))
    standardDeviation = sd(as.numeric(dataSet))
  }
  cvReturn[1] = average
  cvReturn[2] = median
  cvReturn[3] = standardDeviation
  cvReturn[4] = I(standardDeviation / average) * 100
  names(cvReturn) = c("avg", "median", "std", "coeffVar")
  if (is.nan(cvReturn[4]) == T){ #NaN only occurs if its x / 0, where mean is 0 since all cells are 0 valued. Thus set cvReturn as 0 instead.
    cvReturn[4] = 0
  }
  
  return(cvReturn)
}

#Takes a threshold input for determing if the well with the specific bead count clears the threshold. If not, the well must be discounted
#for that sample. This occurs due to the fact that the beads can clump, which may influence the intensity read from the assay, resulting in
#a mischaracterized concentration. 
cleanCounts = function(inputDataSet, threshold){
  cleanData = inputDataSet
  for (i in 1:nrow(inputDataSet)){
    for (j in 3:ncol(inputDataSet)){
      if(as.numeric(inputDataSet[i,j]) >= threshold){
        cleanData[i,j] = as.numeric(inputDataSet[i,j])
      } else {
        cleanData[i,j] = "Remove" #Set cell value as remove 
      }
    }
  }
  return(cleanData)
}

#Takes an input data set and examines for two string characters < and >. These two string characters indicate extremeties in terms of result
#when calculating the relevant concentration from the standard curves. As such, if there is a vast majority of < or > in the column, then 
#the column shall be removed from the data set in order to remove low-variance columns such that lower dimensionality can be achieved. Thus,
#a threshold input is also added, where threshold stands for tolerable percentage of < or > within the column. Anything above the threshold 
#indicates that the column shall be removed
cleanExtremeties = function(inputDataSet, classDataSet, refVariables, threshold, logFileName){
  cat("\n", file = logFileName, append=T)
  cat(paste('#---- Cleaning up variables where number of non-varying values (at max or min of standard curve) exceeds input threshold of ', 
            threshold, '% ----# \n', sep=''), file=logFileName, append=T)
  cat('Referenced class or classes are as follows: \n', file=logFileName, append=T)
  for (a in 1:length(refVariables)){
    cat(paste('->', refVariables[a], ' [With total of ', length(unique(eval(parse(text=paste("classDataSet$", refVariables[a], sep=''))))),
              ' classes] \n', sep=''), file=logFileName, append=T)
  }
  
  removeColumnIndex = matrix(NA, ncol = ncol(inputDataSet), nrow = length(refVariables))
  for (a in 1:length(refVariables)){
    classLabels = unique(eval(parse(text=paste("classDataSet$", refVariables[a], sep=''))))
    cat(paste("\nRemoval based off the current class label: ", refVariables[a], '\n', sep=''), file=logFileName, append=T)
    for (j in 3:ncol(inputDataSet)){
      refLength = length(classLabels) #Retrieves number of classes for current label
      lowcolumnCount = array(0, refLength)
      highcolumnCount = array(0, refLength)
      names(highcolumnCount) = classLabels
      names(lowcolumnCount) = classLabels
        
      for (i in 1:nrow(inputDataSet)){
        if(is.interger0(grep('>', resultDataSet[i,j])) == FALSE){
          currentClass = as.character(eval(parse(text=paste("classDataSet$", refVariables[a], "[i]", sep='')))) #Retrieves class for instance
          highcolumnCount[currentClass] = highcolumnCount[currentClass] + 1
        } else if ((is.interger0(grep('<', resultDataSet[i,j])) == FALSE)){
          currentClass = as.character(eval(parse(text=paste("classDataSet$", refVariables[a], "[i]", sep='')))) #Retrieves class for instance
          lowcolumnCount[currentClass] = lowcolumnCount[currentClass] + 1
        } else {
          lowcolumnCount = lowcolumnCount 
          highcolumnCount = highcolumnCount
        }
      }
      
      obsPerClassArray = array(0, length(classLabels)) #Initializes array to store number of observations per class for proper thresholding
      for (k in 1:length(classLabels)){
        obsPerClassArray[k] = eval(parse(text=paste("sum(classDataSet$", refVariables[a], "== classLabels[k])", sep='')))
      }
      
      #If the requirements where the number of extremeties on one end is greater than the threshold for the whole set, determine if it is 
      #specific to a specific class or if it is random. If random, discard variable.
      if ( any(I(100*lowcolumnCount/obsPerClassArray) > threshold) | any(I(100*highcolumnCount/obsPerClassArray) > threshold) ){
        removeColumnIndex[j-2] = j
        cat("\n", file=logFileName, append=T)
        cat(paste("Will remove the following variable: ", names(resultDataSet)[j], sep=''), file=logFileName, append=T)
        if ( any(I(100*lowcolumnCount/obsPerClassArray) > threshold) & !any(I(100*highcolumnCount/obsPerClassArray) > threshold)){
          cat("\n Reason is because it exceeded threshold for lower bound number of extremities", file=logFileName, append=T)
        } else if ( any(I(100*highcolumnCount/obsPerClassArray) > threshold) & !any(I(100*lowcolumnCount/obsPerClassArray) > threshold))  {
          cat("\n Reason is because it exceeded threshold for higher bound number of extremities", file=logFileName, append=T)
        } else {
          cat("\n Reason is because it exceeded threshold for both lower and higher bound number of extremities", file=logFileName, append=T)
        }
      }
    }
  }
  
  removeColumnIndex = removeColumnIndex[!is.na(removeColumnIndex)]
  cat(paste("\n \nRemoved a total of ", length(removeColumnIndex), " variables", sep=''), file=logFileName, append=T)
  cat("\n#---- Done indexing out variables to be removed ----# \n", file=logFileName, append=T)
  if (length(removeColumnIndex) > 0){
    return(removeColumnIndex)
  } else {
    return(NA)
  }
}

#Finds the cells where NAN is present 
cleanNAN = function(resultDataSet){
  for (i in 3:ncol(resultDataSet)){
    columnCount = 0
    for (j in 1:nrow(resultDataSet)){
      if(is.interger0(grep('N/A', resultDataSet[j,i])) == FALSE){
        resultDataSet[j,i] = NA
      } else if (is.interger0(grep('NA', resultDataSet[j,i], ignore.case=TRUE)) == FALSE){
        resultDataSet[j,i] = NA
      } else if (is.interger0(grep('n/a', resultDataSet[j,i])) == FALSE){
        resultDataSet[j,i] = NA
      } else if (is.interger0(grep('NAN', resultDataSet[j,i], ignore.case=TRUE)) == FALSE){
        resultDataSet[j,i] = NA
      }
    }
  }
  return(resultDataSet)
}

#Function to combine multiple preprocessed datasets and compile it into one dataset with conditions where NA cells/rows/columns 
#are removed. Then it converts all values to numeric automatically as well. It outputs a list, first list is the total merged matrix,
#second object in list is a merged type of file, but does not have string characters with the data tags (e.g. count low), instead it only
#has values. 
mergeSets = function(dataSet1, dataSet2, naCol, logFileName){
  cat("\n", file=logFileName, append=T)
  cat("#---- Merging bead count data set with result data set ----# \n" , file=logFileName, append=T)
  cat("Merging cleaned counts data set with cleaned result data set. Output will be a dataset where low bead counts are logged, and variables 
      that exceed stated threshold of number of extremities per class are removed \n" , file=logFileName, append=T)
  cat("Values that are at limits (i.e. < x or > x), are modified where max limit values take the maximum value detectable, whilst
      minimum limit values are set to 0. \n" , file=logFileName, append=T)
  if (any(is.na(naCol) == F)){ #Checks if there are indices which designate that columns were set to NA previously
    dataSet1 = dataSet1[, -naCol]
    dataSet2 = dataSet2[, -naCol]
  }
  dataSet1_untagged = dataSet1
  for (i in 1:nrow(dataSet1)){
    for (j in 3:ncol(dataSet1)){
      if(is.na(dataSet1[i,j]) == F){
        if(is.interger0(grep('>', dataSet1[i,j])) == FALSE){
          highLimit = strsplit(dataSet1[i,j], '>')
          dataSet1[i,j] = as.numeric(trimws(highLimit[[1]][2]))
        } else if ((is.interger0(grep('<', dataSet1[i,j])) == FALSE)){
          dataSet1[i,j] = 0
        } else {
          dataSet1[i,j] = as.numeric(as.character(dataSet1[i,j]))
        }
      } else {
        dataSet1[i,j] = -100 #Sets NA values as -100, since all relevant values are positive, negative values will be removed subsequently
      }
      dataSet1_untagged[i,j] = dataSet1[i,j]
      
      if (is.interger0(grep('Remove', dataSet2[i,j])) == FALSE){
        dataSet1[i,j] = 'CountLow'
        if (dataSet1_untagged[i,j] == -100){ dataSet1_untagged[i,j] = "CountLow" } 
        else { dataSet1_untagged[i,j] = dataSet1_untagged[i,j] } #Remains the same for the untagged version of the matrix.
      }
    }
  }
  
  #Added for loop to remove leading whitespace in Sample column
  for (i in 1:nrow(dataSet1)){
    dataSet1[i, "Sample"] = trimws(dataSet1[i, "Sample"], which = c("left"))
  }
  for (i in 1:nrow(dataSet1_untagged)){
    dataSet1_untagged[i, "Sample"] = trimws(dataSet1_untagged[i, "Sample"], which = c("left"))
  }
  
  dataSet1 = dataSet1[, !(colSums(dataSet1 < 0) > 0)] #drops columns that have summation less than 0, since NA's were converted to negative.
  dataSet1_untagged =  dataSet1_untagged[, !(colSums(dataSet1_untagged < 0) > 0)] #Same application as above line to the untagged matrix
  cat("\nRemaining variables kept are (excluding location, sample and total-events from the count): \n", file=logFileName, append=T)
  for (i in 3:I(length(names(dataSet1))-1)){
    cat(paste("->", names(dataSet1)[i], "\n", sep=''), file=logFileName, append=T)
  }
  cat(paste("Total number of variables kept are: ", length(names(dataSet1))-3,sep=''), file=logFileName, append=T)
  cat("\n#----Merging Done----# \n", file=logFileName, append=T)
  return(list(dataSet1, dataSet1_untagged))
}

#Function to remove the Background, standards and quality controls.User only has to input the dataset to clean values from. Another input
#is the column index for the user to set the function to refer to for determining standards, background and quality related cells.
cleanNonExperimentVals = function(dataSet, sampleIndex, logFileName){
  backgroundVals = unique(grep('Background', dataSet[[sampleIndex]], ignore.case=T, value=T))
  standardVals = unique(grep('Standard', dataSet[[sampleIndex]], ignore.case=T, value=T))
  qualityVals = unique(grep('Q', dataSet[[sampleIndex]], ignore.case=T, value=T))
  
  
  #Sets background vals cells and rows as NA for removal
  for (i in 1:length(backgroundVals)){
    valIndex = which(backgroundVals[i] == dataSet[[sampleIndex]])
    for (j in 1:length(valIndex)){
      dataSet[valIndex[j], ] = NA #Sets entire row as NA for removal downstream
    }
  }
  
  #Sets standard vals cells and rows as NA for removal
  for (i in 1:length(standardVals)){
    valIndex = which(standardVals[i] == dataSet[[sampleIndex]])
    for (j in 1:length(valIndex)){
      dataSet[valIndex[j], ] = NA #Sets entire row as NA for removal downstream
    }
  }
  
  if (length(qualityVals) > 0){
    #Sets quality control vals cells and rows as NA for removal
    for (i in 1:length(qualityVals)){
      valIndex = which(qualityVals[i] == dataSet[[sampleIndex]])
      for (j in 1:length(valIndex)){
        dataSet[valIndex[j], ] = NA #Sets entire row as NA for removal downstream
      }
    }
    logSentence3 = paste("Removed quality control values from matrix, total of ", length(qualityVals), " unique observation(s) \n", sep='')
  } else {
    logSentence3 = paste("No quality control values detected \n", sep='')
  }
  
  dataSet = dataSet[!is.na(dataSet[[2]]), ] #Removes all NA from dataSet from the sample column
  logSentence1 = paste("Removed background values from matrix, total of ", length(backgroundVals), " unique observation(s) \n", sep='')
  logSentence2 = paste("Removed standard values from matrix, total of ", length(standardVals), " unique observation(s) \n", sep='')
  cat(logSentence1, file=logFileName, append=T)
  cat(logSentence2, file=logFileName, append=T)
  cat(logSentence3, file=logFileName, append=T)
  
  return(dataSet)
}

#Function that examines two datasets, and takes an input string. Using dataset2, also known as the ref matrix, dataset1 will become dataset2
#but will keep its original values only when the string primer occurs for particulars cells in dataset2. Thus, ultimately, dataset1 will
#become dataset2 excluding the cells in dataset2 that have the reference string. 
formMatFromRefMatrix =function(dataSet, refSet, primer){
  for (i in 1:nrow(dataSet)){
    for(j in 1:ncol(dataSet)){
      if (refSet[i,j] == primer){
        dataSet[i,j] = dataSet[i,j] #Remains the same without changing 
      } else {
        dataSet[i,j] = refSet[i,j] #Takes up the reference data set value since 
      }
    }
  }
  return(dataSet)
}

#Function to compress Matrix. Can take in reference indices and exclude Indices as numeric or character names (both have to have the same type)
compressMatrix = function(dataSet, excludeIndices, referenceIndex, logFileName){
  if (is.character(excludeIndices) == T & is.character(referenceIndex) == T){
    referenceIndex = which(names(dataSet) == referenceIndex)
    tempExclude = array(NA, length(excludeIndices))
    for (i in 1:length(excludeIndices)){
      tempExclude[i] = which(names(dataSet) == excludeIndices[i])
    }
    excludeIndices = tempExclude
  }
  logSentence1 = paste("Compressed Class Labels with the following reference Vector: ", names(dataSet)[referenceIndex], "\n", sep='')
  logSentence2 = paste("Compressed Class Labels excluding the following vectors: ", names(dataSet)[excludeIndices], "\n", sep='')
  refVector = dataSet[, referenceIndex]
  refVector = unique(refVector)
  
  dataSetVar = dataSet[, -excludeIndices]
  tempOutput = data.frame(matrix(NA, nrow = length(refVector), ncol = ncol(dataSetVar)))
  names(tempOutput) = names(dataSetVar)
  for (i in 1:length(refVector)){
    subset = dataSetVar[dataSet[,referenceIndex] == refVector[i],]
    subset = subset[1,] #Take first row only, supposedly all rows have equal values across columns 
    tempOutput[i,] = subset
  }
  if (is.na(logFileName) == F){
    cat(logSentence1, file=logFileName, append=T)
    cat(logSentence2, file=logFileName, append=T)
  }
  
  return(tempOutput)
}

#Function to retrieve mode
getmode = function(array){
  unq = unique(array)
  unq[which(tabulate(match(array, unq)) == max(tabulate(match(array, unq))))]
}

#Function for calculating averages from a cleaned data set
calcAverageMat = function(processedDataSet, sampleIndex, removeCols, primer, mode, logFileName){
  subjectVal = unique(processedDataSet[[sampleIndex]])
  tempMat = matrix(NA, nrow = length(subjectVal), ncol = ncol(processedDataSet))
  tempMat = as.data.frame(tempMat)
  colnames(tempMat) = colnames(processedDataSet)
  labelCols = processedDataSet[, c(removeCols)]
  tempMat = tempMat[, -removeCols]
  logSentence1 = paste("Calculating the mode: ", mode, "\n", "Will compress the class matrix as well in the process \n", sep='')
  if (mode == "Mean"){
    for (a in 1:length(subjectVal)){
      currentSubject = subjectVal[a]
      tempMatrix = processedDataSet[processedDataSet[[sampleIndex]] == currentSubject,]
      tempMatrix = tempMatrix[,-removeCols]
      if(nrow(tempMatrix) > 0){
        for (b in 1:ncol(tempMatrix)){
          retrievedMat = calculateCoeffVar(tempMatrix[,b], primer)
          if (is.na(retrievedMat["avg"]) == F){
            tempMat[a, b] = retrievedMat["avg"] #Retrieves the average for that column
          } else {
            tempMat[a,b] = "Low Bead Count Across All Replicates"
          }
        }
      } else {
        print("ERROR IN FINDING AVERAGES/MEDIAN")
      }
    }
  } else if (mode == "Median"){
    for (a in 1:length(subjectVal)){
      currentSubject = subjectVal[a]
      tempMatrix = processedDataSet[processedDataSet[[sampleIndex]] == currentSubject,]
      tempMatrix = tempMatrix[,-removeCols]
      if(nrow(tempMatrix) > 0){
        for (b in 1:ncol(tempMatrix)){
          retrievedMat = calculateCoeffVar(tempMatrix[,b], primer)
          if (is.na(retrievedMat["median"]) == F){
            tempMat[a, b] = retrievedMat["median"] #Retrieves the average for that column
          } else {
            tempMat[a,b] = "Low Bead Count Across All Replicates"
          }
        }
      } else {
        print("ERROR IN FINDING AVERAGES/MEDIAN")
    }
      }
  }
  labelCols = compressMatrix(dataSet = labelCols, excludeIndices = "Location", referenceIndex = "Sample", logFileName)
  tempMat = cbind(labelCols, tempMat) #Created a CV matrix for each individual result
  
  if (is.na(logFileName) == F){
    cat(logSentence1, file=logFileName, append=T)
  }
  return(tempMat)
  
}

#Function that tags an untagged Matrix. All it does is take a matrix, and compresses it, and only denotes if the cells for each set of replicates
#had a string label attached, and will simply return that in the compressed matrix, otherwise it'll return the value "numeric". 
tagUntaggedMat = function(untaggedMat, compressMetric, removeCols, primer){
  tempMat = matrix(NA, nrow=length(unique(untaggedMat[,compressMetric])), ncol=ncol(untaggedMat[,-removeCols]))
  compressedMetricVals = unique(untaggedMat[,compressMetric])
  for (i in 1:nrow(tempMat)){
    subset = untaggedMat[which(untaggedMat[, compressMetric] == compressedMetricVals[i]),-removeCols]
    subsetFinal = array(NA, ncol(subset))
    for (j in 1:ncol(subset)){
      for (k in 1:nrow(subset)){
        if (any(primer == subset[k,j]) == T){
          subsetFinal[j] = subset[k,j] #Keep Value there
          break
        } else {
          subsetFinal[j] = "Optimal"
        }
      }
    }
    
    for (m in 1:length(subsetFinal)){
      tempMat[i,m] = subsetFinal[m]
    }
  }
  return(tempMat)
}

#Returns the magnitude of a vector 
magnitude_vector = function(inputVar){
  return(sqrt(sum(inputVar^2)))
}

#Automatically does z-score standardization 
z_score_standardization = function(inputVar){
  colMeanVals = colMeans(inputVar)
  colSDVals = colSds(as.matrix(inputVar))
  for (j in 1:ncol(inputVar)){
    for (i in 1:nrow(inputVar)){
      inputVar[i,j] = I(inputVar[i,j] - colMeanVals[j]) / colSDVals[j]
    }
  }
  return(inputVar)
}

#---- Debugging Scripts ----

# processedSet = cleanDataOffCVMat(dataSet = resultDataSet, 
#                                  sampleIndex = length(removeColVals), 
#                                  beadCount = beadDataSet,
#                                  removeCols = removeColVals,
#                                  tolerance = 25,  
#                                  beadCountThresh = 20, 
#                                  primerString = "CountLow",
#                                  logFileName = logFile)

#temp = cleanExtremeties(resultDataSet, classesDataSet, cleanRefLabels, 60, logFile) 
#mergedMatrixList = mergeSets(dataSet1 = resultDataSet,
#                            dataSet2 = cleanbeadDataSet, 
#                            naCol = temp, 
#                            logFileName = logFile) #Merges cleaned counts with cleaned results to remove bad counts 
#dgd = tagUntaggedMat(merged_untaggedDataSet, "Sample", c(1,2,3,4,5,6), c("CountLow", "All Low Count, Calc Avg", "Outlier"))
#ded = cleanDataOffCVMat(resultDataSet, sampleIndex = 6, beadDataSet, c(1,2,3,4,5,6), tolerance = 25, beadCountThresh = 20, primerString = "CountLow", logFileName)
#dfd = calcAverageMat(ded, 6, c(1,2,3,4,5,6), c("CountLow", "All Low Count, Calc Avg", "Outlier"), mode="Mean", logFileName)