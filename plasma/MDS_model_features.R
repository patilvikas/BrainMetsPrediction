suppressPackageStartupMessages({
  library(tidyverse)
  library(doParallel)
  library(viridis)
  library(Rtsne)
  library(data.table)
  library(gridExtra)
})

Cluster <-makeCluster(5)
Cluster <- registerDoParallel(Cluster)

Model_out <- "./Models/"
model_name <- "BM"
n_features <- 300
n_times <- 50

MeDIP_RPKM_Data <- readRDS("./data/MeDIPs.rds")

pheno <- fread("./data/SampleInfo.txt")
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

#View a map of the features here
ModelList <- AllIterations_onevEach
ModelList<- do.call(c, unlist(ModelList, recursive=FALSE))
Features <- lapply(ModelList, function(x) x$finalModel$xNames)
Features <- unique(unlist(Features))

#Plot an MDS
MDS <- cmdscale(dist(t(MeDIP_log2CPM[rownames(MeDIP_log2CPM) %in% Features,])),
                k = 2)
MDS <- data.frame(MDS)
MDS <- cbind(MDS, pheno)
annotationBar <- plasma(n_groups)

p1 <- qplot(data = MDS , 
            y = X1, 
            x = X2 , 
            colour = group, 
            size = I(2))+
  geom_point(data = MDS, 
             aes(y = X1, x = X2), 
             size = I(2))+
  theme_bw()+
  scale_colour_manual(values = annotationBar) +
  theme(panel.border = element_blank(), 
  panel.grid.minor = element_blank(), 
  legend.position = "bottom",
  axis.text.x = element_text(colour = "black", size = 11), 
  axis.text.y = element_text(colour = "black", size = 11)) +
  ylab("Dim 1")+
  xlab("Dim 2")
  
pdf(paste0(Model_out, model_name, "_", n_times, "_Iterations_", n_features, "_features.pdf"), 
    height = 4, width = 5)
print(p1)
dev.off()

