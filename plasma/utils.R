suppressPackageStartupMessages({
  library(tidyverse)
  library(reshape2)
  library(caret)
  library(glmnet)
  library(doParallel)
  library(pROC)
  library(viridis)
  library(data.table)
  library(gridExtra)
})

SplitFunction <- function(Mat, Classes)
{
  
  require(dplyr)
  require(caret)
  require(glmnet)
  
  ##Split into training and test sets
  
  df <- data.frame(ID = colnames(Mat),
                   Classes = Classes)
  samples <- createDataPartition(df$Classes, 
                                 p = 0.8,
                                 times = 50)
  
  return(list(df = df, samples = samples))
}


PredFunction <- function(ModelList, TestData, Indices, classes.df)
{ 
  TrainPheno <- classes.df[Indices,]
  TestData <- TestData[,!colnames(TestData) %in% TrainPheno$ID]
  TestPheno <- classes.df%>%filter(!ID %in% TrainPheno$ID)
  
  Predictions.list <- list()
  OutputNames <- names(ModelList)
  
  for(i in 1:length(ModelList))
  {
    
    Features <- ModelList[[i]]$Model$finalModel$xNames
    TestDataNew <- TestData[match(Features,rownames(TestData)), ]
    Model <- ModelList[[i]]$Model
    Prediction.classProbs <- predict(Model, 
                                     newdata = t(TestDataNew), 
                                     type = "prob") %>%
      data.frame
    
    Prediction.classProbs$ActualClass <- TestPheno$Classes
    Prediction.classProbs$PredictedClass <- predict(Model, 
                                                    newdata = t(TestDataNew),
                                                    type = "raw")
    
    Predictions.list[[i]] <- Prediction.classProbs
    message(i)
  }
  
  names(Predictions.list) <- OutputNames
  return(Predictions.list)
  
}
