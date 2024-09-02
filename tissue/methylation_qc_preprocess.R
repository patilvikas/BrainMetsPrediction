library(minfi)
library(maxprobes)
library(tidyverse)
library(data.table)
library(RColorBrewer)
library(limma)

###################################################
### adjust path and read idata files
###################################################
idatPath <- "data/IDATS/"

#analysisPath
analysisPath <- "data/Results/"
ifelse(!dir.exists(analysisPath), dir.create(analysisPath), FALSE)
setwd(analysisPath)

project <- "BrainMets"
###################################################
### Read sample info
###################################################
targets <- fread(paste0(idatPath, "SampleSheet.csv"), 
                 skip = 7)
targets <- as.data.frame(targets)
targets$IDAT <- paste0(targets$Sentrix_ID, "_",  targets$Sentrix_Position)
targets$Basename <- file.path(idatPath, targets$IDAT)

###################################################
### Read idat files
###################################################
RGset <- read.metharray.exp(targets = targets,
                            force=TRUE)
sampleNames(RGset) <- targets$Sample_Name

###################################################
######  calculate the detection p-values
###################################################
detP <- detectionP(RGset)
detPV <- colMeans(detP)
names(detPV) <- targets$Sample_Name
detPV <- as.data.frame(detPV)
colnames(detPV) <- "Average_detection_pval"
detPV <- detPV %>%
  tibble::rownames_to_column("Sample_Name")

pal <- brewer.pal(8,"Dark2")
pdf(paste0(project, "_QC_Allsample_detection_p-value.pdf"), 
    width = 6, height = 5)
par(mfrow=c(1,1))
barplot(colMeans(detP), 
        col=pal[factor(targets$Sentrix_ID)], 
        las=2,
        cex.names=0.8,
        ylab="Mean detection p-values")
abline(h=0.05,col="red")
dev.off()

###################################################
### predict material type FFPE or Frozen
###################################################
getFFPE <- function(RGset,thr=2000){
  control_probes <- getProbeInfo(RGset, type = "Control")
  green <- getGreen(RGset)
  idx <- which(rownames(green) %in% control_probes[control_probes$Type=="RESTORATION","Address"])
  ifelse(green[idx,]>thr,"FFPE","Frozen")
}

materialType <- getFFPE(RGset)
materialType <- as.data.frame(materialType)
materialType <- materialType %>%
  tibble::rownames_to_column("Sample_Name")

###################################################
### normalization 
###################################################
Mset.noob <- preprocessNoob(RGset)

#### plot density plot
pdf(paste0(project, "_density_plot.pdf"),
    width = 5.5, height = 4)
for(i in 1:dim(Mset.noob)[2]){
  densityPlot(Mset.noob[,i],
              main="Normalized", 
              legend=FALSE,
              pal = "black")
  legend("top", 
         legend = levels(factor(targets$Sample_Name[i])),
         text.col= "black")
}
dev.off()

###################################################
### MDS plots to look at largest sources of variation
###################################################
pdf(paste0(project, "_MDS_plot_before-filtering.pdf"), 
    width = 8, height = 5)
plotMDS(getM(Mset.noob), 
        top=1000,
        gene.selection="common",
        col=pal[factor(targets$Sentrix_ID)])
dev.off()

###################################################
### get ratioset
###################################################
Rset <- ratioConvert(Mset.noob, 
                     what = "both",
                     keepCN=T)
GRset <- mapToGenome(Rset)
# add SNP info
RSetSNP <- addSnpInfo(GRset)

###################################################
### get Beta values
###################################################
beta_value <- getBeta(RSetSNP)
beta_value <- as.data.frame(beta_value)

###################################################
### LUMP_Purity
###################################################
LUMP_44_CpG_site<- c("cg00240653","cg00450164","cg00880290","cg00933696","cg01138020",
                     "cg02026204","cg02053964","cg02167021","cg02997560","cg03431741","cg03436397",
                     "cg03841065","cg04915566","cg05199874","cg05305434","cg05769344","cg05798125",
                     "cg07002058","cg07598052","cg07641284","cg08854008","cg09302355","cg09606470",
                     "cg10511890","cg10559416","cg13030790","cg13912307","cg14076977","cg14913777",
                     "cg17518965","cg19466818","cg20170223","cg20695297","cg21164509","cg21376733",
                     "cg22331159","cg23114964","cg23553480","cg24796554","cg25384897","cg25574765",
                     "cg26427109","cg26842802","cg27215100")

LUMP_CpG <- beta_value[rownames(beta_value) %in% LUMP_44_CpG_site, ]
LUMP <- colMeans(LUMP_CpG)/0.85
LUMP <- as.data.frame(LUMP)
LUMP <- LUMP %>%
  tibble::rownames_to_column("Sample_Name")

###################################################
### compile QC results
###################################################
FinalQcResult <- list(targets, detPV, materialType, LUMP) %>% 
  reduce(inner_join)
FinalQcResult <- FinalQcResult %>%
  dplyr::select(-c(Basename))

write.table(FinalQcResult, 
            file = paste0(project, "_Meth_QC_results.txt"), 
            quote =FALSE, sep = "\t", eol = "\n", dec = ".",  
            row.names = FALSE, col.names = TRUE)

saveRDS(beta_value, 
        paste0(project, "_allSamples_allProbes_beatvalues.rds"))
