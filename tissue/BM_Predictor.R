library(data.table)
library(tidyverse)
library(caret)

outpath <- "./data/Results/Predictor_scores/"
ifelse(!dir.exists(outpath), dir.create(outpath), FALSE)

load("./data/model/GBMobjs.RData")

beta_value <- readRDS("./data/Results/TestingCohort_allProbes_beatvalues.rds")
beta_value <- as.data.frame(beta_value)

### check all probes are present in testing cohort
common_probes <- intersect(rownames(beta_value), predictorNames)
if(length(predictorNames) != length(common_probes)){
  stop("features are not matching!!")
}
####

newDF1 <- t(beta_value[rownames(beta_value)%in%predictorNames,])

predicted_probs <- predict(objGBM, newDF1, type='prob')
rownames(predicted_probs) <- rownames(newDF1)
colnames(predicted_probs) <- paste0(colnames(predicted_probs), "_prob")
predicted_probs <- predicted_probs %>%
  tibble::rownames_to_column("Sample_Name")

predicted_probs <- predicted_probs[order(predicted_probs$yes_prob, decreasing = T), ]

####
write.table(predicted_probs,
            paste0(outpath, "Tumor_meth_BM_predicted_probabilities.txt"),
            sep = "\t", eol = "\n", quote = F,
            row.names = F, col.names = T)
