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

source("utils.R")

Cluster <-makeCluster(5)
Cluster <- registerDoParallel(Cluster)

Model_out <- "./Models/"
model_name <- "BM"
n_features <- 300
n_times <- 50

MeDIP_RPKM_Data <- readRDS("./data/MeDIPs.rds")

pheno <- fread("./data/Testing_SampleInfo.txt")
pheno$Sample <- pheno$SampleID
pheno <- pheno %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames(var="Sample") %>% 
  as.data.frame

MeDIP_RPKM_Data <- MeDIP_RPKM_Data[, rownames(pheno)]
MeDIP_log2CPM <- log2(MeDIP_RPKM_Data * 0.3 + 1e-6)
groups <- as.vector(unique(pheno$group))
n_groups <- length(unique(pheno$group))

#load model
load(paste0(Model_out, model_name, "_", n_times, "_Iterations_", n_features, "_features.RData"))

#load split information
load(paste0(Model_out, model_name, "_Splits.RData"))

Classes.df <- Splits$df
TestPerformance_list <- list()
for(i in 1:n_times)
{
  TestPerformance_list[[i]] <- PredFunction(ModelList = AllIterations_onevEach[[i]],
                                            TestData = MeDIP_log2CPM, 
                                            Indices = Splits$samples[[i]], 
                                            classes.df = Classes.df)
}

save(TestPerformance_list, 
     file = paste0(Model_out, model_name, "_TestPerformance_", n_times, "_Iterations_", n_features, "_features.RData") )


