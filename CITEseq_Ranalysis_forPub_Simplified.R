library(tibble)
library(reticulate)
library(sctransform)
library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(glmGamPoi)
library(patchwork)
library(Azimuth)
library(SeuratWrappers)
library(celldex)
library(SingleR)
library(viridis)
library(SingleCellExperiment)
library(scater)
library(scran)
library(patchwork)
library(MAST)
library(msigdbr)
library(enrichR)
library(fgsea)
library(tidyr) 
library(cowplot)
library(edgeR)
library(Matrix)
library(S4Vectors)
library(pheatmap)
library(apeglm)
library(DESeq2)
library(data.table)
library(purrr)
library(scWidgets)
library(EnhancedVolcano)
library(Seurat)
library(SeuratObject)
library(SeuratData)
library(SeuratDisk)
library(SeuratWrappers)


mypalette <- brewer.pal(12,"Paired")
mypalette2 <- brewer.pal(12,"Set3")
mypalette3 <- brewer.pal(11,"RdBu")

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}




## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


##########################

setwd("./Seurat5")

##Check cell cycle gene expression
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#### NOTE THAT ADT COCKTAIL FOR EXP 2, 3, 4, 5 used the 220+ antibody panel while EXP 6, 7, 8 used the updated 160+ antibody panel. Unfortunately, the 220 antibody cocktail was no longer available through Biolegend.

############.  EXP 6 ##################
# read raw data using the Seurat function "Read10X" 
cells = Seurat::Read10X_h5(filename = "MDS_Exp6_CellBender_filtered_seurat.h5",use.names = T)

#Identify singlets first using HTO information
GEX.data <- cells$"Gene Expression" 
ADT.HTO.data <- cells$"Antibody Capture"
ADT.data <- ADT.HTO.data[-grep("HTO",rownames(ADT.HTO.data)),] #split ADT
HTO.data <- ADT.HTO.data[c("HTO1","HTO2","HTO3","HTO4","HTO5","HTO6"),] #split ADT
pbmc <- CreateSeuratObject(counts = GEX.data, project = "MDS_exp6", min.cells = 3, min.features = 200) 
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
joint.bcs <- intersect(colnames(pbmc), colnames(ADT.data))
ADT.data <- ADT.data[, joint.bcs]
HTO.data <- HTO.data[, joint.bcs]
pbmc[["ADT"]] <- CreateAssayObject(counts=ADT.data)
pbmc[["HTO"]] <- CreateAssayObject(counts=HTO.data)
pbmc <- NormalizeData(pbmc, assay = "HTO", normalization.method = "CLR")
pbmc <- HTODemux(pbmc, assay = "HTO", positive.quantile = 0.99)
table(pbmc$HTO_classification.global)
#Doublet Negative  Singlet 
#6769      485    16105 
s_exp6 <- subset(pbmc, subset = nFeature_RNA >=200 & nCount_ADT >500  & HTO_classification.global == "Singlet" & percent.mt <= 10) #10% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8599307/#:~:text=In%20summary%2C%20we%20reported%20a,seq%20QC%20of%20human%20samples.
table(s_exp6$HTO_classification.global)
#Doublet Singlet 
#5837   14029 
prots_new <- rownames(pbmc[["ADT"]]) # for ADT integration below


############.  EXP 7 ##################
# read raw data using the Seurat function "Read10X" 
cells = Seurat::Read10X_h5(filename = "MDS_Exp7_CellBender_filtered_seurat.h5",use.names = T)

#Identify singlets first using HTO information
GEX.data <- cells$"Gene Expression" 
ADT.HTO.data <- cells$"Antibody Capture"
ADT.data <- ADT.HTO.data[-grep("HTO",rownames(ADT.HTO.data)),] #split ADT
HTO.data <- ADT.HTO.data[c("HTO7","HTO8","HTO9","HTO10","HTO12","HTO13"),] #split ADT
pbmc <- CreateSeuratObject(counts = GEX.data, project = "MDS_exp7", min.cells = 3, min.features = 200)
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
joint.bcs <- intersect(colnames(pbmc), colnames(ADT.data))
ADT.data <- ADT.data[, joint.bcs]
HTO.data <- HTO.data[, joint.bcs]
pbmc[["ADT"]] <- CreateAssayObject(counts=ADT.data)
pbmc[["HTO"]] <- CreateAssayObject(counts=HTO.data)
pbmc <- NormalizeData(pbmc, assay = "HTO", normalization.method = "CLR")
pbmc <- HTODemux(pbmc, assay = "HTO", positive.quantile = 0.99)
table(pbmc$HTO_classification.global)
# Doublet Negative  Singlet 
#5359      537    11232 

s_exp7 <- subset(pbmc, subset = nFeature_RNA >=200 & nCount_ADT >500  & HTO_classification.global == "Singlet" & percent.mt <= 10) #10% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8599307/#:~:text=In%20summary%2C%20we%20reported%20a,seq%20QC%20of%20human%20samples.
table(s_exp7$HTO_classification.global)
#Doublet Singlet 
#4057    8248


############.  EXP 7 Brep2 (to get more reads for samples) ##################
# read raw data using the Seurat function "Read10X" 
cells = Seurat::Read10X_h5(filename = "MDS_Exp7b_CellBender_filtered_seurat.h5",use.names = T)

#Identify singlets first using HTO information
GEX.data <- cells$"Gene Expression" 
ADT.HTO.data <- cells$"Antibody Capture"
ADT.data <- ADT.HTO.data[-grep("HTO",rownames(ADT.HTO.data)),] #split ADT
HTO.data <- ADT.HTO.data[c("HTO7","HTO8","HTO9","HTO10","HTO12","HTO13"),] #split ADT
pbmc <- CreateSeuratObject(counts = GEX.data, project = "MDS_exp7b", min.cells = 3, min.features = 200)
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
joint.bcs <- intersect(colnames(pbmc), colnames(ADT.data))
ADT.data <- ADT.data[, joint.bcs]
HTO.data <- HTO.data[, joint.bcs]
pbmc[["ADT"]] <- CreateAssayObject(counts=ADT.data)
pbmc[["HTO"]] <- CreateAssayObject(counts=HTO.data)
pbmc <- NormalizeData(pbmc, assay = "HTO", normalization.method = "CLR")
pbmc <- HTODemux(pbmc, assay = "HTO", positive.quantile = 0.99)
table(pbmc$HTO_classification.global)
#Doublet Negative  Singlet 
#7502     3324    10652 
s_exp7b <- subset(pbmc, subset = nFeature_RNA >=200 & nCount_ADT >500  & HTO_classification.global == "Singlet" & percent.mt <= 10) #10% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8599307/#:~:text=In%20summary%2C%20we%20reported%20a,seq%20QC%20of%20human%20samples.
table(s_exp7b$HTO_classification.global)
#Doublet Singlet 
#6154    8641 



############.  EXP 8 ##################
# read raw data using the Seurat function "Read10X" 
cells = Seurat::Read10X_h5(filename = "MDS_Exp8_CellBender_filtered_seurat.h5",use.names = T)

#Identify singlets first using HTO information
GEX.data <- cells$"Gene Expression" 
ADT.HTO.data <- cells$"Antibody Capture"
ADT.data <- ADT.HTO.data[-grep("HTO",rownames(ADT.HTO.data)),] #split ADT
HTO.data <- ADT.HTO.data[c("HTO1","HTO2","HTO3","HTO4","HTO5","HTO6"),] #split ADT
pbmc <- CreateSeuratObject(counts = GEX.data, project = "MDS_exp8", min.cells = 3, min.features = 200)
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
joint.bcs <- intersect(colnames(pbmc), colnames(ADT.data))
ADT.data <- ADT.data[, joint.bcs]
HTO.data <- HTO.data[, joint.bcs]
pbmc[["ADT"]] <- CreateAssayObject(counts=ADT.data)
pbmc[["HTO"]] <- CreateAssayObject(counts=HTO.data)
pbmc <- NormalizeData(pbmc, assay = "HTO", normalization.method = "CLR")
pbmc <- HTODemux(pbmc, assay = "HTO", positive.quantile = 0.99)
table(pbmc$HTO_classification.global)
#Doublet Negative  Singlet 
#2823      227    11927
s_exp8 <- subset(pbmc, subset = nFeature_RNA >=200 & nCount_ADT >500  & HTO_classification.global== "Singlet" & percent.mt <= 10)
table(s_exp8$HTO_classification.global)
#Singlet 
#11066 



############  EXP 2  ##################
# read raw data using the Seurat function "Read10X" 
cells = Seurat::Read10X_h5(filename = "MDS_Exp2_CellBender_filtered_seurat.h5",use.names = T)

#Identify singlets first using HTO information
GEX.data <- cells$"Gene Expression" 
ADT.HTO.data <- cells$"Antibody Capture"
ADT.data <- ADT.HTO.data[-grep("HTO",rownames(ADT.HTO.data)),] #split ADT
HTO.data <- ADT.HTO.data[c("HTO6","HTO7","HTO8","HTO9","HTO10"),] #split ADT
pbmc <- CreateSeuratObject(counts = GEX.data, project = "MDS_exp2", min.cells = 3, min.features = 200)
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
joint.bcs <- intersect(colnames(pbmc), colnames(ADT.data))
ADT.data <- ADT.data[, joint.bcs]
HTO.data <- HTO.data[, joint.bcs]
pbmc[["ADT"]] <- CreateAssayObject(counts=ADT.data)
pbmc[["HTO"]] <- CreateAssayObject(counts=HTO.data)
pbmc <- NormalizeData(pbmc, assay = "HTO", normalization.method = "CLR")
pbmc <- HTODemux(pbmc, assay = "HTO", positive.quantile = 0.99)
table(pbmc$HTO_classification.global)
#Doublet Negative  Singlet 
#2473     2705     9455
s_exp2 <- subset(pbmc, subset = nFeature_RNA >=200 & nCount_ADT >500  & HTO_classification.global== "Singlet" & percent.mt <= 10)
table(s_exp2$HTO_classification.global)
#Singlet 
#8373


############  EXP 3  ##################
# read raw data using the Seurat function "Read10X" 
cells = Seurat::Read10X_h5(filename = "MDS_Exp3_CellBender_filtered_seurat.h5",use.names = T)

#Identify singlets first using HTO information
GEX.data <- cells$"Gene Expression" 
ADT.HTO.data <- cells$"Antibody Capture"
ADT.data <- ADT.HTO.data[-grep("HTO",rownames(ADT.HTO.data)),] #split ADT
HTO.data <- ADT.HTO.data[c("HTO1","HTO2","HTO3","HTO4","HTO5"),] #split ADT
pbmc <- CreateSeuratObject(counts = GEX.data, project = "MDS_exp3", min.cells = 3, min.features = 200)
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
joint.bcs <- intersect(colnames(pbmc), colnames(ADT.data))
ADT.data <- ADT.data[, joint.bcs]
HTO.data <- HTO.data[, joint.bcs]
pbmc[["ADT"]] <- CreateAssayObject(counts=ADT.data)
pbmc[["HTO"]] <- CreateAssayObject(counts=HTO.data)
pbmc <- NormalizeData(pbmc, assay = "HTO", normalization.method = "CLR")
pbmc <- HTODemux(pbmc, assay = "HTO", positive.quantile = 0.99)
table(pbmc$HTO_classification.global)

#Doublet Negative  Singlet 
#9877     8320    19724
s_exp3 <- subset(pbmc, subset = nFeature_RNA >=200 & nCount_ADT >500  & HTO_classification.global== "Singlet" & percent.mt <= 10)
table(s_exp3$HTO_classification.global)
#Singlet 
# 15648 



############  EXP 4  ##################
# read raw data using the Seurat function "Read10X" 
cells = Seurat::Read10X_h5(filename = "MDS_Exp4_CellBender_filtered_seurat.h5",use.names = T)

#Identify singlets first using HTO information
GEX.data <- cells$"Gene Expression" 
ADT.HTO.data <- cells$"Antibody Capture"
ADT.data <- ADT.HTO.data[-grep("HTO",rownames(ADT.HTO.data)),] #split ADT
HTO.data <- ADT.HTO.data[c("HTO6","HTO7","HTO8","HTO9","HTO10"),] #split ADT
pbmc <- CreateSeuratObject(counts = GEX.data, project = "MDS_exp4", min.cells = 3, min.features = 200) #change min cells to account for issues with cell cycle scoring
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
joint.bcs <- intersect(colnames(pbmc), colnames(ADT.data))
ADT.data <- ADT.data[, joint.bcs]
HTO.data <- HTO.data[, joint.bcs]
pbmc[["ADT"]] <- CreateAssayObject(counts=ADT.data)
pbmc[["HTO"]] <- CreateAssayObject(counts=HTO.data)
pbmc <- NormalizeData(pbmc, assay = "HTO", normalization.method = "CLR")
pbmc <- HTODemux(pbmc, assay = "HTO", positive.quantile = 0.99)
table(pbmc$HTO_classification.global)
#Doublet Negative  Singlet 
#3910     4741    14841 
s_exp4 <- subset(pbmc, subset = nFeature_RNA >=200 & nCount_ADT >500  & HTO_classification.global== "Singlet" & percent.mt <= 10)
table(s_exp4$HTO_classification.global)
#Singlet 
#12820  



############  EXP 5  ##################
# read raw data using the Seurat function "Read10X" 
cells = Seurat::Read10X_h5(filename = "MDS_Exp5_CellBender_filtered_seurat.h5",use.names = T)

#Identify singlets first using HTO information
GEX.data <- cells$"Gene Expression" 
ADT.HTO.data <- cells$"Antibody Capture"
ADT.data <- ADT.HTO.data[-grep("HTO",rownames(ADT.HTO.data)),] #split ADT
HTO.data <- ADT.HTO.data[c("HTO1","HTO2","HTO3","HTO4","HTO5"),] #split ADT
pbmc <- CreateSeuratObject(counts = GEX.data, project = "MDS_exp5", min.cells = 3, min.features = 200)
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
joint.bcs <- intersect(colnames(pbmc), colnames(ADT.data))
ADT.data <- ADT.data[, joint.bcs]
HTO.data <- HTO.data[, joint.bcs]
pbmc[["ADT"]] <- CreateAssayObject(counts=ADT.data)
pbmc[["HTO"]] <- CreateAssayObject(counts=HTO.data)
pbmc <- NormalizeData(pbmc, assay = "HTO", normalization.method = "CLR")
pbmc <- HTODemux(pbmc, assay = "HTO", positive.quantile = 0.99)
table(pbmc$HTO_classification.global)
#Doublet Negative  Singlet 
#4010     2571    12123 
s_exp5 <- subset(pbmc, subset = nFeature_RNA >=200 & nCount_ADT >500  & HTO_classification.global== "Singlet" & percent.mt <= 10)
table(s_exp5$HTO_classification.global)
#Singlet 
#9872



remove(ADT.HTO.data)
remove(GEX.data)
remove(HTO.data)
remove(pbmc)


######## INTEGRATE DIFFERENT EXPERIMENTS INTO ONE DATA SET AND PERFORM SC TRANSFORM ###############
BM_merged <- merge(s_exp2, y = c(s_exp3,s_exp4,s_exp5,s_exp6,s_exp7,s_exp7b,s_exp8), add.cell.ids = c("Exp2", "Exp3","Exp4","Exp5", "Exp6","Exp7","Exp7b","Exp8"), project = "BMmerged")
BM_merged$HTO_sampleinfo <- paste(BM_merged$orig.ident, BM_merged$HTO_classification,sep="-")

filteredmeta <- data.frame(rownames(BM_merged@meta.data),BM_merged@meta.data[,"HTO_sampleinfo"])
write.table(filteredmeta, "BMmerged_metadata_HTOclassification_mito10.txt", sep = "\t", col.names = F, row.names = F, quote = F) #Add clinical metadata to this file
metaclinical <- read.table("BMmerged_metadata_HTOclassification_wClinical.txt", sep = "\t", header = F, stringsAsFactors = F)
rownames(metaclinical) <- metaclinical$V1
metaclinical <- metaclinical[,c(3:16)]
colnames(metaclinical) <- c("SampleType","Time","Time2","Patient","Institution","Age","Gender","Race","IPSS",	"TotalCyclesGiven",	"BestResponse",	"moOS",	"moPFS","Survival")

BM_merged<-AddMetaData(object = BM_merged, metadata = metaclinical)

SaveSeuratRds(BM_merged, "BM_merged_AllExperiments.rds")

#Use SCtransform to integrate#
DefaultAssay(BM_merged) <-'RNA'
options(future.globals.maxSize= 891289600) #address global size isue
BM.list <- SplitObject(BM_merged, split.by = "orig.ident")# cant do it on every sample due to matrix error (too many vectors)
BM.list.sct <- lapply(X = BM.list, FUN = SCTransform, method = "glmGamPoi", vars.to.regress = c("percent.mt"), vst.flavor = "v2") ## use SCTransform v2 to allow DEG analysis
BM.list.sct <- lapply(X = BM.list.sct, FUN = CellCycleScoring, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
for (i in 1:length(BM.list.sct)) {BM.list.sct[[i]]$CC.Difference <- BM.list.sct[[i]]$S.Score - BM.list.sct[[i]]$G2M.Score}
for (i in 1:length(BM.list.sct)) {DefaultAssay(BM.list.sct[[i]]) <-'RNA'}
BM.list.sct <- lapply(X = BM.list.sct, FUN = SCTransform, method = "glmGamPoi", vars.to.regress = c("CC.Difference","percent.mt"), vst.flavor = "v2") ## regress out cell cycle difference and percent.mt

features <- SelectIntegrationFeatures(object.list = BM.list.sct, nfeatures = 3000)
features <- features[-c(grep("MALAT",features),grep("^MT-",features))] #remove Mito and MALAT1 genes as anchors for integration
BM.list.sct <- PrepSCTIntegration(object.list = BM.list.sct, anchor.features = features)
BM.list.sct <- lapply(X = BM.list.sct, FUN = RunPCA, features = features)

immune.anchors <- FindIntegrationAnchors(object.list = BM.list.sct, normalization.method = "SCT",reference = c(1,2,3,4,5,6,8), #reference the experiment with most cells (ie. exp 2,3 (BM), exp5 (PBMC), and exp6,7,8 (BM))
                                         anchor.features = features, dims = 1:40, reduction = "rpca", k.anchor = 20)
all_features <- lapply(BM.list.sct , row.names) %>% Reduce(intersect, .) 

BM.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

BM.combined.sct <- RunPCA(BM.combined.sct, verbose = FALSE)
BM.combined.sct <- RunUMAP(BM.combined.sct, reduction = "pca", dims = 1:40)


#### Normalize genes from RNA assay ####
DefaultAssay(BM.combined.sct) = "RNA"
BM.combined.sct<-NormalizeData(BM.combined.sct)
all.genes <- rownames(BM.combined.sct)
BM.combined.sct <- ScaleData(BM.combined.sct, features = all.genes)
DefaultAssay(BM.combined.sct) <- "integrated"



remove(BM.list)
remove(BM.list.sct)
remove(immune.anchors)
remove(s_exp2)
remove(s_exp3)
remove(s_exp4)
remove(s_exp5)
remove(s_exp6)
remove(s_exp7)
remove(s_exp7b)
remove(s_exp8)




######################     Integrate ADT data #########
DefaultAssay(BM_merged) <- "ADT"
ADT.list <- SplitObject(BM_merged, split.by = "orig.ident")
ADT.list <- lapply(X = ADT.list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "CLR", margin = 2)
  FindVariableFeatures(x,  selection.method = "vst", verbose = FALSE)
})
features <- prots_new[-grep("ctrl",prots_new)] # remove IGG control as feature
ADT.list <- lapply(X = ADT.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})


ADT.anchors <- FindIntegrationAnchors(object.list = ADT.list, anchor.features = features, reduction = "rpca",k.anchor = 20, dims = 1:20)
ADT.combined <- IntegrateData(anchorset = ADT.anchors,dims = 1:20) 

DefaultAssay(ADT.combined) <- "integrated"

ADT.combined <- ScaleData(ADT.combined, verbose = FALSE)
ADT.combined <- RunPCA(ADT.combined, verbose = FALSE)
ADT.combined <- RunUMAP(ADT.combined, reduction = "pca", dims = 1:20)
ADT.combined <- FindNeighbors(ADT.combined, reduction = "pca", dims = 1:20)
ADT.combined <- FindClusters(ADT.combined, resolution = 0.5)


DimPlot(ADT.combined, reduction = "umap", group.by = "orig.ident", label = F,
        repel = TRUE)
FeaturePlot(ADT.combined, features = c("CD34.1","CD117-c-kit","CD19.1","CD20","CD45RO","CD45RA","CD4.1","CD8","CD314-NKG2D","CD335-NKp46" ,"CD14.1","CD16","CD11b","CD11c","CD105","HLA-DR"),
            reduction = 'umap', min.cutoff = "q50",max.cutoff = "q95", ncol = 4)& scale_colour_gradientn(colours = c("light gray", brewer.pal(n = 9, name = "YlOrRd")))

remove(ADT.list)
remove(ADT.anchors)



#######LOAD SCtransformed + ADT INTEGRATED DATA#####################
BM.combined.sct[["integrated.adt"]]<-ADT.combined[["integrated"]]
BM.combined.sct[["adt.pca"]]<-ADT.combined[["pca"]]
BM.combined.sct[["ADT2"]]<-ADT.combined[["ADT"]]

ADT.combined <- NULL



###run Seurat integration and WNN 
BM.combined.sct = FindMultiModalNeighbors(  
  BM.combined.sct, reduction.list = list("pca", "adt.pca"), 
  dims.list = list(1:40, 1:20),
  modality.weight.name = "RNA.weight",
  verbose = FALSE
)


BM.combined.sct  <- RunUMAP(BM.combined.sct, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

no

BM.combined.sct<- FindClusters(BM.combined.sct, graph.name = "wsnn", 
                               algorithm = 3, 
                               resolution = 3, 
                               verbose = FALSE)





###########Visualize clusters #########

pp1 <- DimPlot(BM.combined.sct, reduction = "wnn.umap", label = F,group.by = "wsnn_res.3",raster=FALSE) + NoLegend() 
LabelClusters(pp1, id = "wsnn_res.3",  fontface = "bold", color = "black", size =4)

f<-FeaturePlot(BM.combined.sct, features = "nCount_RNA",reduction = 'wnn.umap', min.cutoff = "200", max.cutoff = "8000",raster=FALSE)& scale_colour_gradientn(colours = c("light gray", brewer.pal(n = 9, name = "Greens")[4:9]))
f1<-FeaturePlot(BM.combined.sct, features = "nFeature_RNA",reduction = 'wnn.umap', min.cutoff = "200", max.cutoff = "4000",raster=FALSE)& scale_colour_gradientn(colours = c("light gray", brewer.pal(n = 9, name = "Greens")[4:9]))
f2<-FeaturePlot(BM.combined.sct, features = "nCount_ADT",reduction = 'wnn.umap', min.cutoff = "500", max.cutoff = "q95",raster=FALSE)& scale_colour_gradientn(colours = c("light gray", brewer.pal(n = 9, name = "Reds")[4:9]))
f3<-FeaturePlot(BM.combined.sct, features = "nFeature_ADT",reduction = 'wnn.umap', min.cutoff = "200", max.cutoff = "q99",raster=FALSE)& scale_colour_gradientn(colours = c("light gray", brewer.pal(n = 9, name = "Reds")[4:9]))
f4<-FeaturePlot(BM.combined.sct, features = "nCount_HTO",reduction = 'wnn.umap', min.cutoff = "100", max.cutoff = "5000",raster=FALSE)& scale_colour_gradientn(colours = c("light gray", brewer.pal(n = 9, name = "Reds")[4:9]))
f5<-FeaturePlot(BM.combined.sct, features = "percent.mt",reduction = 'wnn.umap', min.cutoff = "0", max.cutoff = "15",raster=FALSE)& scale_colour_gradientn(colours = c("light gray", brewer.pal(n = 9, name = "Purples")[4:9]))

patchwork::wrap_plots(f,f1,f2,f5, ncol = 2) 



#####################FIND ALL MARKERS FOR CLUSTERS ##########################
BM.combined.sct1<-PrepSCTFindMarkers(BM.combined.sct, assay = "SCT", verbose = T) # need to run this before performing FindMarkers or DEG calling

Idents(BM.combined.sct1) = "wsnn_res.3"
BM.combined.sct.markers <- FindAllMarkers(BM.combined.sct1, assay = 'SCT',only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#ADT differential Expresison

Idents(BM.combined.sct1) = "wsnn_res.3"
BM.combined.sct.ADTmarkers <- FindAllMarkers(BM.combined.sct1, assay = 'ADT2',only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)



############ RUN Azimuth with Bone Marrow Reference #####
available_data <- AvailableData()

DefaultAssay(BM.combined.sct) <- "SCT"
BM.combined.sct[["SCT"]] <- split(BM.combined.sct[["SCT"]], f = BM.combined.sct$orig.ident ) #not sure why but need to re-split and joinlayers the SCT for azimuth to work
BM.combined.sct[["SCT"]] <- JoinLayers(BM.combined.sct[["SCT"]])

bm <- RunAzimuth(BM.combined.sct, reference = "bonemarrowref", assay = "SCT") #issues with RefDR so downloaded the reference data from azimuth
p1<-DimPlot(bm, reduction = "wnn.umap", group.by = "predicted.celltype.l2", label = TRUE, repel = TRUE, label.size = 3, raster=F, order=T) 

p4<-DimPlot(bm, reduction = "wnn.umap", group.by =  "predicted.celltype.l2", cols = c('CD14 Mono'='#ffc5d9',
                                                                                      'CD16 Mono'='#c2f2d0',
                                                                                      'CD4 Memory'='#349c44',
                                                                                      'CD4 Naive'='#e0a899',
                                                                                      'CD56 bright NK'='#f3ea5f',
                                                                                      'CD8 Effector_1'='#ed5555',
                                                                                      'CD8 Effector_2'='#b4c468',
                                                                                      'CD8 Memory'='#43e8d8',
                                                                                      'CD8 Memory_1'='#b8a7ea',
                                                                                      'CD8 Memory_2'='#b7e0dc',
                                                                                      'CD8 Naive'='blue',
                                                                                      'cDC1'='#71c7ec',
                                                                                      'cDC2'='#ff3f3f',
                                                                                      'ydT'='#c04df9',
                                                                                      'GMP'='#1b85b8',
                                                                                      'CLP'='#00ff00',
                                                                                      'HSC'='#800000',
                                                                                      'LMPP'='#ffcc5c',
                                                                                      'Early Eryth'='#66b2b2',
                                                                                      'Late Eryth'='#006666',
                                                                                      'BaEoMa'="blue",
                                                                                      'Platelet'='#189ad3',
                                                                                      'MAIT'='#6b3e26',
                                                                                      'Memory B'='#ae5a41',
                                                                                      'transitional B'='#fe00f6',
                                                                                      'Naive B'='#c3cb71',
                                                                                      'NK'='#baffc9',
                                                                                      'pDC'='purple',
                                                                                      'Plasma'='#bae1ff',
                                                                                      'Prog_B 1'='#ff71ce',
                                                                                      'Prog_B 2'='#01cdfe',
                                                                                      'Prog_DC'='#ffcb85',
                                                                                      'Prog Mk'='#000000',
                                                                                      'EMP'='#f8d3c5',
                                                                                      'Treg'='#fbad26'),
            
            label = F, repel = TRUE, label.size = 3, raster=F, order=F)
LabelClusters(p4, id = "predicted.celltype.l2",  fontface = "bold", color = "black", size =6)









bm2<- RunAzimuth(BM.combined.sct, reference = "pbmcref", assay = "SCT")
ppp4<-DimPlot(bm2, reduction = "wnn.umap", group.by =  "predicted.celltype.l2", cols = c('ASDC'='#ffc5d9',
                                                                                         'B intermediate'='#c2f2d0',
                                                                                         'B memory'='#349c44',
                                                                                         'B naive'='#e0a899',
                                                                                         'CD14 Mono'='#f3ea5f',
                                                                                         'CD16 Mono'='#ed5555',
                                                                                         'CD4 CTL'='#fe00f6',
                                                                                         'CD4 Naive'='#0bff01',
                                                                                         'CD4 Proliferating'='#b7e0dc',
                                                                                         'CD4 TCM'='#bdd1a0',
                                                                                         'CD4 TEM'='#ff3f3f',
                                                                                         'CD8 Naive'='#c04df9',
                                                                                         'CD8 Proliferating'='#1b85b8',
                                                                                         'CD8 TCM'='#5a5255',
                                                                                         'CD8 TEM'='#bb88dd',
                                                                                         'cDC1'='#a2a2d0',
                                                                                         'cDC2'='#ae5a41',
                                                                                         'dnT'='#ff91af',
                                                                                         'Eryth'='#bd7657',
                                                                                         'ydT'='#586b41',
                                                                                         'HSPC'='#bae1ff',
                                                                                         'ILC'='#ff71ce',
                                                                                         'MAIT'='#01cdfe',
                                                                                         'NK'='#ffcb85',
                                                                                         'NK Proliferating'='#000000',
                                                                                         'NK_CD56bright'='#345a54',
                                                                                         'pDC'='#fbad26',
                                                                                         'Plasmablast'='#00008b',
                                                                                         'Platelet'='#0e1111',
                                                                                         'Treg'='#0057e7'),
              
              label = F, repel = F, label.size = 3, raster=F,order=F)
LabelClusters(ppp4, id = "predicted.celltype.l2",  fontface = "bold", color = "black", size =6)





############## Make DotPLOT visualizing the markers using well-established cell type-specific markers for cluster identification ################
###   Broad markers      ####
broadfeatures_gex <- c("CD34","PDGFRA","ENG","TFRC","ITGA2B","IL3RA","FCER1A","CEACAM8","NCAM1","TRGV9","KLRB1","CD3E","CD3G","CD8A","CD4","CD14","FCGR3A","ITGAM","ITGAX","CD1C","CD19","MS4A1")
broadfeatures_adt <- c("CD34.1","CD140a-PDGFRa","CD105","CD71","CD41","CD123","FceRIa","CD66b","CD56","TCR-Vy9","TCR-Vd2","CD161", "TCR-Va7.2","CD8","CD4.1","CD45RO","CD45RA","CD14.1","CD16","CD11b","CD11c","CD1c","CD19.1","CD20")


Idents(BM.combined.sct) = "wsnn_res.3"
p<-DotPlot(BM.combined.sct, cluster.idents = T, group.by =  "manualAnno_broad", features = broadfeatures_gex, assay = "SCT",scale.min=0,scale.max=60)+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+ scale_colour_gradient2(low = "white", mid = "#71c7ec", high = "#107dac",midpoint=1)+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p1<-DotPlot(BM.combined.sct,cluster.idents = T, group.by =  "manualAnno_broad", features = broadfeatures_adt, assay = "ADT2",scale.min=0,scale.max=100)+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+ scale_colour_gradient2(low = "white", mid = "#fceee9", high = "#ff0000",midpoint=0.5)+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
patchwork::wrap_plots(p,p1, ncol = 2)

broadfeatures_gex2 <- c("CD34","PDGFRA","ENG","TFRC","HBB","CA1","JCHAIN","STMN1","HBG2","AHSP","TUBB","AREG","ITGA2B","IL3RA","FCER1A","CEACAM8","NCAM1","TRGV9","KLRB1","CD3E","CD3G","CD8A","CD4","CD14","FCGR3A","ITGAM","ITGAX","CD1C","CD19","MS4A1","IGLV2-14","IGKV3-15","IGKV1-5")
DotPlot(BM.combined.sct,cluster.idents=T, features = broadfeatures_gex2, assay = "SCT",scale.min=0,scale.max=60)+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+ scale_colour_gradient2(low = "white", mid = "#71c7ec", high = "#107dac",midpoint=1)+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

broadfeatures_gex3 <- c("CD34","PDGFRA","ENG","TFRC","HBB","CA1","JCHAIN","STMN1","HBG2","AHSP","AREG","ITGA2B","IL3RA","S100A8","S100A9","FCER1A","CEACAM8","NCAM1","NKG7","TRGV9","KLRB1","CD3E","CD3G","CD8A","CD4","CD14","FCGR3A","ITGAM","ITGAX","AIF1","C1QA","C1QB","CD1C","CD19","MS4A1","IGLV2-14","IGKV3-15","IGKV1-5")
d<-DotPlot(BM.combined.sct,cluster.idents=T, features = broadfeatures_gex3, assay = "SCT",scale.min=0,scale.max=60)+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+ scale_colour_gradient2(low = "white", mid = "#71c7ec", high = "#107dac",midpoint=1)+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

DotPlot(BM.combined.sct,cluster.idents = T, features = c("CD4.1","CD8","CD45RO","CD45RA","CD197-CCR7","CD62L","CD28.1","CD127-IL-7Ra","CD44.1","KLRG1-MAFA"), assay = "ADT2",scale.min=0,scale.max=100)+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+ scale_colour_gradient2(low = "white", mid = "#fceee9", high = "#ff0000",midpoint=0.5)+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))




################################ MANUAL VERIFICATION of CLUSTER ANNOTATION after CAREFUL REVIEW ##############################
manualAnno <- read.table("BMmerged_WNN3_manual_cluster_annotation.txt", header = T, sep = "\t", stringsAsFactors = F)
BM.combined.sct@meta.data$wsnn_res.3 <- droplevels(BM.combined.sct@meta.data$wsnn_res.3)
BM.combined.sct@meta.data$manualAnno_broad <- NA
for (i in unique(BM.combined.sct@meta.data$wsnn_res.3)) {
  BM.combined.sct@meta.data[BM.combined.sct@meta.data$wsnn_res.3 == i,]$manualAnno_broad <-  manualAnno[manualAnno$Cluster ==i,]$Broad
}
BM.combined.sct@meta.data$manualAnno_broad2 <- NA
for (i in unique(BM.combined.sct@meta.data$wsnn_res.3)) {
  BM.combined.sct@meta.data[BM.combined.sct@meta.data$wsnn_res.3 == i,]$manualAnno_broad2 <-  manualAnno[manualAnno$Cluster ==i,]$Broad2
}




Idents(BM.combined.sct) = "wsnn_res.3"
BM.combined.sct<- subset(x = BM.combined.sct, idents = c(18,55,61,64,70,91), invert = TRUE) ## Remove low quality or artifact clusters deteremined by marker gene or ADT 

BM.combined.sct$manualAnno_broad <- factor(BM.combined.sct$manualAnno_broad, levels = rev(c("CD34","EPC","Megakaryocyte","Neutrophil","Myeloid","T & NK","B cell")))

pp<-DimPlot(BM.combined.sct, reduction = "wnn.umap", group.by =  "manualAnno_broad",cols = c('B cell'='#6565bf',
                                                                                             'Myeloid'='#ffe700',
                                                                                             "Megakaryocyte"='grey',
                                                                                             'EPC'='#dbac98',
                                                                                             'CD34'='#ff0000',
                                                                                             'T & NK'='#71c7ec',
                                                                                             'Neutrophil'='#0e1111'),
            
            label = F, repel = F, label.size = 3, raster=F,order=T)+ NoLegend() 
LabelClusters(pp, id = "manualAnno_broad",  fontface = "bold", color = "black", size =6)+ NoLegend() 


broadfeatures_gex <- c("CD34","ENG","TFRC","ITGA2B","FCER1A","CEACAM8","CD14","ITGAM","ITGAX","CD1C","FCGR3A","NCAM1","CD4","CD8A","CD19","MS4A1")
broadfeatures_adt <- c("CD34.1","CD105","CD71","CD41","FceRIa","CD66b","CD14.1","CD11b","CD11c","CD1c","CD16","CD56","CD4.1","CD8","CD19.1","CD20")

p<-DotPlot(BM.combined.sct, cluster.idents = F, group.by =  "manualAnno_broad", features = broadfeatures_gex, assay = "SCT",scale.min=0,scale.max=60)+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+ scale_colour_gradient2(low = "white", mid = "#71c7ec", high = "#107dac",midpoint=1)+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p1<-DotPlot(BM.combined.sct,cluster.idents = F, group.by =  "manualAnno_broad", features = broadfeatures_adt, assay = "ADT2",scale.min=0,scale.max=100)+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+ scale_colour_gradient2(low = "white", mid = "#fceee9", high = "#ff0000",midpoint=0.5)+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
patchwork::wrap_plots(p,p1, ncol = 2)




pp1<-DimPlot(BM.combined.sct, reduction = "wnn.umap", group.by =  "manualAnno_broad2",cols = c('B cell'='#6565bf',
                                                                                               'Plasmablast'='#b4a7d6',
                                                                                               'Myeloid'='#ffe700',
                                                                                               'CD4'='#99cc33',
                                                                                               'CD8'='#02a9f7',
                                                                                               'EPC'='#dbac98',
                                                                                               'gdT'='#e7b416',
                                                                                               'CD34'='#ff0000',
                                                                                               'MAIT'='#7df9ff',
                                                                                               'NK'='#f6b26b',
                                                                                               'Neutrophil'='#0e1111',
                                                                                               'Megakaryocyte'='#00467F'),
             
             label = F, repel = F, label.size = 3, raster=F,order=T)
LabelClusters(pp1, id = "manualAnno_broad2",  fontface = "bold", color = "black", size =6)










################ SUBSET BROAD ANNOTATIONS TO RECLUSTER AND IDENTIFY SPECIFIC ANNOTAIONS
Idents(BM.combined.sct) = "manualAnno_broad"
BM.combined.sct_myelo <- subset(x = BM.combined.sct, idents = c("CD34","EPC","Myeloid"))
BM.combined.sct_TNK <- subset(x = BM.combined.sct, idents = c("T & NK"))
BM.combined.sct_Bcell <- subset(x = BM.combined.sct, idents = c("B cell"))
BM.combined.sct_Others <- subset(x = BM.combined.sct, idents = c("CD34","EPC","Myeloid","T & NK","B cell"), invert = TRUE)


#Annotate cells in "Others (Neutro/mega)" group
BM.combined.sct_Others$manualAnno_specific  <- BM.combined.sct_Others$manualAnno_broad2
BM.combined.sct_Others$manualAnno_specific2  <- BM.combined.sct_Others$manualAnno_broad2


##### need to rerun SCTransform and re-integrate data for myelo data
DefaultAssay(BM.combined.sct_myelo) <- "RNA"
BM.combined.sct_myelo.list <- SplitObject(object = BM.combined.sct_myelo , split.by = "orig.ident")
BM.combined.sct_myelo.list <- lapply(X = BM.combined.sct_myelo.list, FUN = SCTransform, method = "glmGamPoi", vars.to.regress = c("CC.Difference","percent.mt"), vst.flavor = "v2") ## regress out cell cycle difference and percent.mt

features <- SelectIntegrationFeatures(object.list = BM.combined.sct_myelo.list, nfeatures = 3000)
features <- features[-c(grep("MALAT",features),grep("^MT-",features))] #remove Mito and MALAT1 genes as anchors for integration
BM.combined.sct_myelo.list <- PrepSCTIntegration(object.list = BM.combined.sct_myelo.list, anchor.features = features)
BM.combined.sct_myelo.list <- lapply(X = BM.combined.sct_myelo.list, FUN = RunPCA, features = features)
immune.anchors <- FindIntegrationAnchors(object.list = BM.combined.sct_myelo.list, normalization.method = "SCT",reference = c(1,2,3,4,5,6,8), #reference the experiment with most cells (ie. exp 2,3 (BM), exp5 (PBMC), and exp6,7,8 (BM))
                                         anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
all_features <- lapply(BM.combined.sct_myelo.list , row.names) %>% Reduce(intersect, .) 
myelo.integrated<- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

myelo.integrated<- RunPCA(myelo.integrated, verbose = FALSE)
myelo.integrated<- RunUMAP(myelo.integrated, reduction = "pca", dims = 1:30)

#re-integrate ADT
DefaultAssay(BM.combined.sct_myelo) <- "ADT"
ADT.list <- SplitObject(BM.combined.sct_myelo, split.by = "orig.ident")
ADT.list <- lapply(X = ADT.list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "CLR", margin = 2)
  x <-FindVariableFeatures(x,  selection.method = "vst", verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = ADT.list)
features <- features[-grep("ctrl",features)] # remove IGG control as feature
ADT.list <- lapply(X = ADT.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
ADT.anchors <- FindIntegrationAnchors(object.list = ADT.list, anchor.features = features, reduction = "rpca",k.anchor = 20, dims = 1:20)
ADT.combined <- IntegrateData(anchorset = ADT.anchors,dims = 1:20) 

DefaultAssay(ADT.combined) <- "integrated"

ADT.combined <- ScaleData(ADT.combined, verbose = FALSE)
ADT.combined <- RunPCA(ADT.combined, verbose = FALSE)

myelo.integrated[["integrated.adt"]]<-ADT.combined[["integrated"]]
myelo.integrated[["adt.pca"]]<-ADT.combined[["pca"]]
myelo.integrated[["ADT2"]]<-ADT.combined[["ADT"]]

myelo.integrated = FindMultiModalNeighbors(  
  myelo.integrated, reduction.list = list("pca", "adt.pca"), 
  dims.list = list(1:30, 1:20),
  modality.weight.name = "RNA.weight",
  verbose = FALSE
)
myelo.integrated <- RunUMAP(myelo.integrated , nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

myelo.integrated <- FindClusters(myelo.integrated, graph.name = "wsnn", 
                                 algorithm = 3, 
                                 resolution = 3, 
                                 verbose = FALSE)


pp <- DimPlot(myelo.integrated, reduction = "wnn.umap", label = F,group.by = "wsnn_res.3",raster=FALSE) + NoLegend() 
LabelClusters(pp, id = "wsnn_res.3",  fontface = "bold", color = "black", size =6)



#annotate clusters after careful review of marker gene expression
manualAnno <- read.table("BMmerged_WNN3_manual_cluster_annotation_CD34_Myeloid_EPC.txt", header = T, sep = "\t", stringsAsFactors = F)
myelo.integrated@meta.data$wsnn_res.3 <- droplevels(myelo.integrated@meta.data$wsnn_res.3)
myelo.integrated@meta.data$manualAnno_specific <- NA
for (i in unique(myelo.integrated@meta.data$wsnn_res.3)) {
  myelo.integrated@meta.data[myelo.integrated@meta.data$wsnn_res.3 == i,]$manualAnno_specific <-  manualAnno[manualAnno$Cluster ==i,]$Specific
}
myelo.integrated@meta.data$manualAnno_specific2 <- NA
for (i in unique(myelo.integrated@meta.data$wsnn_res.3)) {
  myelo.integrated@meta.data[myelo.integrated@meta.data$wsnn_res.3 == i,]$manualAnno_specific2 <-  manualAnno[manualAnno$Cluster ==i,]$Specific2
}

Idents(myelo.integrated) <- "manualAnno_specific2"

myelo.integrated$manualAnno_specific2 <- factor(myelo.integrated$manualAnno_specific2, levels = rev(c("EPC","MEP","HSC","LMPP","GMP","Classical Mono","Intermediate Mono","Non-classical Mono",'mDC','pDC','Neutrophil',"Stromal",'Differentiating Stem Cell')))

DimPlot(myelo.integrated, reduction = "wnn.umap", group.by =  "manualAnno_specific2",cols = c("HSC"='#ff0000',
        "MEP"='#a70000',
        "LMPP"='#bca9bd',
        "GMP"='pink',
        'Differentiating Stem Cell'='#ffdfba',
        "Stromal"='#7FFFCC',
        "EPC"='#dbac98',
        "Non-classical Mono"='#876127',
        "Classical Mono"='#d0b783',
        "Intermediate Mono"='yellow',
        "mDC"='#c27ba0',
        "pDC"='#326ada','Neutrophil'='darkgrey'), label = F, repel = F, label.size = 3, raster=F,order=T)+ NoLegend() 

p<-DotPlot(myelo.integrated, cluster.idents = F,features = c("ENG","TFRC", "HPGDS","GATA2","CD34","KIT","CRHBP","AVP","NPR3","CD38","NPW","PRTN3","ELANE","AZU1" ,"CD4","CD14","FCGR3A","ITGAM","ITGAX","CD1C","CLEC4C","CEACAM8","CXCL12","CD3E"), assay = "SCT",scale.min=0,scale.max=60)+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+ scale_colour_gradient2(low = "white", mid = "#F5F8FA", high = "#107dac")+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p1<-DotPlot(myelo.integrated,cluster.idents = F, features = c("CD105","CD71","CD34.1","CD117-c-kit","CD4.1","CD14.1","CD16","CD11b","CD11c","CD1c","CD303-BDCA2","HLA-DR","CD66b","CD3"), assay = "ADT2",scale.min=0,scale.max=100)+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+ scale_colour_gradient2(low = "white", mid = "#F8F1F4", high = "#BF3F3F")+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
multiplot(p,p1)




########################### need to rerun SCTransform and re-integrate data for TNKB data ##################
set.seed(123)
DefaultAssay(BM.combined.sct_TNK) <- "RNA"
BM.combined.sct_TNK.list <- SplitObject(object = BM.combined.sct_TNK , split.by = "orig.ident")
BM.combined.sct_TNK.list <- lapply(X = BM.combined.sct_TNK.list, FUN = SCTransform, method = "glmGamPoi", vars.to.regress = c("CC.Difference","percent.mt"), vst.flavor = "v2") ## regress out cell cycle difference and percent.mt

features <- SelectIntegrationFeatures(object.list = BM.combined.sct_TNK.list, nfeatures = 3000)
features <- features[-c(grep("MALAT",features),grep("^MT-",features))] #remove Mito and MALAT1 genes as anchors for integration
BM.combined.sct_TNK.list <- PrepSCTIntegration(object.list = BM.combined.sct_TNK.list, anchor.features = features)
BM.combined.sct_TNK.list <- lapply(X = BM.combined.sct_TNK.list, FUN = RunPCA, features = features)
immune.anchors <- FindIntegrationAnchors(object.list = BM.combined.sct_TNK.list, normalization.method = "SCT",reference = c(1,2,3,4,5,6,8),  
                                         anchor.features = features, dims = 1:30, reduction = "rpca")
all_features <- lapply(BM.combined.sct_TNK.list , row.names) %>% Reduce(intersect, .) 
TNK.integrated<- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

TNK.integrated<- RunPCA(TNK.integrated, verbose = FALSE)
TNK.integrated<- RunUMAP(TNK.integrated, reduction = "pca", dims = 1:30)

#re-integrate ADT
DefaultAssay(BM.combined.sct_TNK) <- "ADT"
ADT.list <- SplitObject(BM.combined.sct_TNK, split.by = "orig.ident")
ADT.list <- lapply(X = ADT.list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "CLR", margin = 2)
  x <-FindVariableFeatures(x,  selection.method = "vst", verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = ADT.list)
features <- features[-grep("ctrl",features)] # remove IGG control as feature
ADT.list <- lapply(X = ADT.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
ADT.anchors <- FindIntegrationAnchors(object.list = ADT.list, anchor.features = features, reduction = "rpca",k.anchor = 20, dims = 1:20)
ADT.combined <- IntegrateData(anchorset = ADT.anchors,dims = 1:20) 

DefaultAssay(ADT.combined) <- "integrated"

ADT.combined <- ScaleData(ADT.combined, verbose = FALSE)
ADT.combined <- RunPCA(ADT.combined, verbose = FALSE)

TNK.integrated[["integrated.adt"]]<-ADT.combined[["integrated"]]
TNK.integrated[["adt.pca"]]<-ADT.combined[["pca"]]
TNK.integrated[["ADT2"]]<-ADT.combined[["ADT"]]

TNK.integrated = FindMultiModalNeighbors(  
  TNK.integrated, reduction.list = list("pca", "adt.pca"), 
  dims.list = list(1:30, 1:20),
  modality.weight.name = "RNA.weight",
  verbose = FALSE
)
TNK.integrated <- RunUMAP(TNK.integrated , nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

TNK.integrated <- FindClusters(TNK.integrated, graph.name = "wsnn", 
                               algorithm = 3, 
                               resolution = 5, 
                               verbose = FALSE) # need higher resolution to separate out CD4 & vd1 gdT from CD8 cells in high mito percent clusters


pp <- DimPlot(TNK.integrated, reduction = "wnn.umap", label = F,group.by = "wsnn_res.5",raster=FALSE) + NoLegend() 
LabelClusters(pp, id = "wsnn_res.5",  fontface = "bold", color = "black", size =6)

###annotate clusters after careful review of marker gene expression
manualAnno <- read.table("BMmerged_WNN5_manual_cluster_annotation_TNK.txt", header = T, sep = "\t", stringsAsFactors = F)
TNK.integrated@meta.data$wsnn_res.5 <- droplevels(TNK.integrated@meta.data$wsnn_res.5)
TNK.integrated@meta.data$manualAnno_specific <- NA
for (i in unique(TNK.integrated@meta.data$wsnn_res.5)) {
  TNK.integrated@meta.data[TNK.integrated@meta.data$wsnn_res.5 == i,]$manualAnno_specific <-  manualAnno[manualAnno$Cluster ==i,]$Specific
}
TNK.integrated@meta.data$manualAnno_specific2 <- NA
for (i in unique(TNK.integrated@meta.data$wsnn_res.5)) {
  TNK.integrated@meta.data[TNK.integrated@meta.data$wsnn_res.5 == i,]$manualAnno_specific2 <-  manualAnno[manualAnno$Cluster ==i,]$Specific2
}


TNK.integrated$manualAnno_specific2 <- factor(TNK.integrated$manualAnno_specific2 , levels = rev(c("CD4 naive","CD4 TCM","CD4 TEM","CD4 CTL","CD4 Treg","CD8 naive","CD8 TCM","CD8 TEM","CD8 TEMRA","MAIT","Vd1 gdT","Vd2 gdT","CD56dim NK","CD56bright NK",'Differentiating Stem Cell')))
Idents(TNK.integrated) <- "manualAnno_specific2"

DimPlot(TNK.integrated, reduction = "wnn.umap", group.by =  "manualAnno_specific2", cols = c("CD4 TCM"="#99cc33",
                                                                                             "CD4 naive"="#339900",
                                                                                             "CD4 Treg"="#ff00c1",
                                                                                             "CD4 CTL"="#8cff32",
                                                                                             "CD4 TEM"="#D6FE83",
                                                                                             "CD8 naive"='#bae1ff',
                                                                                             "CD8 TCM"='#b3cde0',
                                                                                             "CD8 TEM"='#005b96',
                                                                                             "CD8 TEMRA"='#5bc0de',
                                                                                             'MAIT'='#7df9ff',
                                                                                             "Vd1 gdT"='#ffd600',
                                                                                             "Vd2 gdT"='#e7b416',
                                                                                             "CD56dim NK"='#ff9966',
                                                                                             "CD56bright NK"='#ff4f00','Differentiating Stem Cell'='#ffdfba'),label = F, repel = F, label.size = 3, raster=F,order=T)+NoLegend()


p<-DotPlot(TNK.integrated, cluster.idents = F,features = c("CD3E","CD4","CD8A","CCR7","SELL","CD28","FOXP3","CTLA4","TRGV9","TRDV1","TRDV2","KLRB1","FCGR3A","NCAM1","CD34"), assay = "SCT",scale.min=0,scale.max=60)+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+ scale_colour_gradient2(low = "white", mid = "#F5F8FA", high = "#107dac")+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p1<-DotPlot(TNK.integrated,cluster.idents = F, features = c("CD3","CD4.1","CD8","CD45RA","CD45RO","CD62L","CD28.1","CD25","TCRa-b","TCR-Vy9","TCR-Vd2","CD161","CD16","CD56","CD34.1"), assay = "ADT2",scale.min=0,scale.max=100)+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+ scale_colour_gradient2(low = "white", mid = "#F8F1F4", high = "#BF3F3F")+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

multiplot(p,p1)





################# PERFORM B CELL RECLUSTERING ##############
DefaultAssay(BM.combined.sct_Bcell) <- "RNA"
BM.combined.sct_Bcell.list <- SplitObject(object = BM.combined.sct_Bcell , split.by = "orig.ident")
BM.combined.sct_Bcell.list <- lapply(X = BM.combined.sct_Bcell.list, FUN = SCTransform, method = "glmGamPoi", vars.to.regress = c("CC.Difference","percent.mt"), vst.flavor = "v2") ## regress out cell cycle difference and percent.mt

features <- SelectIntegrationFeatures(object.list = BM.combined.sct_Bcell.list, nfeatures = 3000)
features <- features[-c(grep("MALAT",features),grep("^MT-",features))] #remove Mito and MALAT1 genes as anchors for integration
BM.combined.sct_Bcell.list <- PrepSCTIntegration(object.list = BM.combined.sct_Bcell.list, anchor.features = features)
BM.combined.sct_Bcell.list <- lapply(X = BM.combined.sct_Bcell.list, FUN = RunPCA, features = features)
immune.anchors <- FindIntegrationAnchors(object.list = BM.combined.sct_Bcell.list, normalization.method = "SCT",reference = c(1,2,3,4,5,6,8),  
                                         anchor.features = features, dims = 1:20, reduction = "rpca")
all_features <- lapply(BM.combined.sct_Bcell.list , row.names) %>% Reduce(intersect, .) 
Bcell.integrated<- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

Bcell.integrated<- RunPCA(Bcell.integrated, verbose = FALSE)
Bcell.integrated<- RunUMAP(Bcell.integrated, reduction = "pca", dims = 1:20)

#re-integrate ADT
DefaultAssay(BM.combined.sct_Bcell) <- "ADT"
ADT.list <- SplitObject(BM.combined.sct_Bcell, split.by = "orig.ident")
ADT.list <- lapply(X = ADT.list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "CLR", margin = 2)
  x <-FindVariableFeatures(x,  selection.method = "vst", verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = ADT.list)
features <- features[-grep("ctrl",features)] # remove IGG control as feature
ADT.list <- lapply(X = ADT.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
ADT.anchors <- FindIntegrationAnchors(object.list = ADT.list, anchor.features = features, reduction = "rpca",k.anchor = 20, dims = 1:10)
ADT.combined <- IntegrateData(anchorset = ADT.anchors,dims = 1:10) 

DefaultAssay(ADT.combined) <- "integrated"

ADT.combined <- ScaleData(ADT.combined, verbose = FALSE)
ADT.combined <- RunPCA(ADT.combined, verbose = FALSE)

Bcell.integrated[["integrated.adt"]]<-ADT.combined[["integrated"]]
Bcell.integrated[["adt.pca"]]<-ADT.combined[["pca"]]
Bcell.integrated[["ADT2"]]<-ADT.combined[["ADT"]]

Bcell.integrated = FindMultiModalNeighbors(  
  Bcell.integrated, reduction.list = list("pca", "adt.pca"), 
  dims.list = list(1:20, 1:10),
  modality.weight.name = "RNA.weight",
  verbose = FALSE
)
Bcell.integrated <- RunUMAP(Bcell.integrated , nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

Bcell.integrated <- FindClusters(Bcell.integrated, graph.name = "wsnn", 
                                 algorithm = 3, 
                                 resolution = 2, 
                                 verbose = FALSE) # need higher resolution to separate out CD4 from CD8 cells in clusters


pp <- DimPlot(Bcell.integrated, reduction = "wnn.umap", label = F,group.by = "wsnn_res.2",raster=FALSE) + NoLegend() 
LabelClusters(pp, id = "wsnn_res.2",  fontface = "bold", color = "black", size =6)


#annotate clusters after careful review of marker gene expression
manualAnno <- read.table("BMmerged_WNN2_manual_cluster_annotation_Bcell.txt", header = T, sep = "\t", stringsAsFactors = F)
Bcell.integrated@meta.data$wsnn_res.2 <- droplevels(Bcell.integrated@meta.data$wsnn_res.2)
Bcell.integrated@meta.data$manualAnno_specific <- NA
for (i in unique(Bcell.integrated@meta.data$wsnn_res.2)) {
  Bcell.integrated@meta.data[Bcell.integrated@meta.data$wsnn_res.2 == i,]$manualAnno_specific <-  manualAnno[manualAnno$Cluster ==i,]$Specific
}
Bcell.integrated@meta.data$manualAnno_specific2 <- NA
for (i in unique(Bcell.integrated@meta.data$wsnn_res.2)) {
  Bcell.integrated@meta.data[Bcell.integrated@meta.data$wsnn_res.2 == i,]$manualAnno_specific2 <-  manualAnno[manualAnno$Cluster ==i,]$Specific2
}

Bcell.integrated$manualAnno_specific2 <- factor(Bcell.integrated$manualAnno_specific2, levels = rev(c("Naive B","Memory B","Plasmablast",'Differentiating Stem Cell')))
  
Idents(Bcell.integrated) <- "manualAnno_specific2"
DimPlot(Bcell.integrated, reduction = "wnn.umap", group.by =  "manualAnno_specific2", cols = c("Naive B"='#efbbff',
                                                                                                 "Memory B"='#6565bf',
                                                                                                 'Plasmablast'='#660066','Differentiating Stem Cell'='#ffdfba'),label = F, repel = F, label.size = 3, raster=F,order=T)+NoLegend()


p<-DotPlot(Bcell.integrated, cluster.idents = F,features = c("CD19","MS4A1","CD24","IGHM","IGHD","CD27","CD38","CD34"), assay = "SCT",scale.min=0,scale.max=60)+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+ scale_colour_gradient2(low = "white", mid = "#F5F8FA", high = "#107dac")+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p1<-DotPlot(Bcell.integrated,cluster.idents = F, features = c("CD19.1","CD20","CD24.1","IgM","IgD","CD27.1","CD38.1","CD34.1"), assay = "ADT2",scale.min=0,scale.max=100)+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+ scale_colour_gradient2(low = "white", mid = "#F8F1F4", high = "#BF3F3F")+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

multiplot(p,p1)




################################ MERGE ALL SEURAT Metadata  TOGETHER WITH SPECIFIC ANNOTATIONS #########################
meta1 <- myelo.integrated@meta.data[,c("manualAnno_specific","manualAnno_specific2")]
meta2 <- TNK.integrated@meta.data[,c("manualAnno_specific","manualAnno_specific2")]
meta3 <- Bcell.integrated@meta.data[,c("manualAnno_specific","manualAnno_specific2")]
meta4 <- BM.combined.sct_Others@meta.data[,c("manualAnno_specific","manualAnno_specific2")]
meta_specific <- rbind(meta1,meta2,meta3,meta4)

meta_specific$manualAnno_specific2 <- as.character(meta_specific$manualAnno_specific2 )
meta_specific[meta_specific$manualAnno_specific2 == "CD4 naive",]$manualAnno_specific2 <- "CD4 Naive"
meta_specific[meta_specific$manualAnno_specific2 == "CD8 naive",]$manualAnno_specific2 <- "CD8 Naive"

BM.combined.sct<-AddMetaData(object = BM.combined.sct, metadata = meta_specific)


Idents(BM.combined.sct) <- "manualAnno_specific"
pqq<-DimPlot(BM.combined.sct, reduction = "wnn.umap", group.by =  "manualAnno_specific", cols = c('B cell'='#6565bf',
                                                                                                  'Monocyte'='#ffe700',
                                                                                                  'Differentiating Stem Cell'='#ffdfba',
                                                                                                  'CD4'='#99cc33',
                                                                                                  'CD8'='#02a9f7',
                                                                                                  'DC'='#c27ba0',
                                                                                                  "Stromal"='#7FFFCC',
                                                                                                  "Megakaryocyte"='lightgrey',
                                                                                                  'EPC'='#dbac98',
                                                                                                  'gdT'='#e7b416',
                                                                                                  'CD34'='#ff0000',
                                                                                                  'MAIT'='#7df9ff',
                                                                                                  'Plasmablast'='#660066',
                                                                                                  'NK'='#ff7400',
                                                                                                  'Neutrophil'='darkgrey'),
             label = F, repel = F, label.size = 3, raster=F,order=T,ncol=1)
LabelClusters(pqq, id = "manualAnno_specific",  fontface = "bold", color = "black", size =4)



pq<-DimPlot(BM.combined.sct, reduction = "wnn.umap",group.by =  "manualAnno_specific2", cols = c("HSC"='#ff0000',
                                                                                                 "MEP"='#a70000',
                                                                                                 "LMPP"='#bca9bd',
                                                                                                 "GMP"='pink',
                                                                                                 'Differentiating Stem Cell'='#ffdfba',
                                                                                                 "Megakaryocyte"='lightgrey',
                                                                                                 "Stromal"='#7FFFCC',
                                                                                                 "EPC"='#dbac98',
                                                                                                 "Non-classical Mono"='#876127',
                                                                                                 "Classical Mono"='#d0b783',
                                                                                                 "Intermediate Mono"='yellow',
                                                                                                 "mDC"='#c27ba0',
                                                                                                 "pDC"='#326ada',
                                                                                                 "Naive B"='#efbbff',
                                                                                                 "Memory B"='#6565bf',
                                                                                                 'Plasmablast'='#660066',
                                                                                                 "CD4 TCM"="#99cc33",
                                                                                                 "CD4 Naive"="#339900",
                                                                                                 "CD4 Treg"="#ff00c1",
                                                                                                 "CD4 CTL"="#8cff32",
                                                                                                 "CD4 TEM"="#D6FE83",
                                                                                                 "CD8 Naive"='#bae1ff',
                                                                                                 "CD8 TCM"='#b3cde0',
                                                                                                 "CD8 TEM"='#005b96',
                                                                                                 "CD8 TEMRA"='#5bc0de',
                                                                                                 'MAIT'='#7df9ff',
                                                                                                 "Vd1 gdT"='#ffd600',
                                                                                                 "Vd2 gdT"='#e7b416',
                                                                                                 "CD56dim NK"='#ff9966',
                                                                                                 "CD56bright NK"='#ff4f00',
                                                                                                 'Neutrophil'='darkgrey'
),label = F, repel = F, label.size = 3, raster=F,order=F)+ NoLegend()
LabelClusters(pq, id = "manualAnno_specific2",  fontface = "bold", color = "black", size =4)+ theme(legend.position = 'none')
pq+ theme(legend.position = 'none')



########Generate MARKER GENE EXPRESSION FIGURE FOR PAPER ######################################
Idents(BM.combined.sct) = "manualAnno_broad"
levels(BM.combined.sct) <- c("CD34","EPC","Myeloid","T & NK", "B cell","Neutrophil","Megakaryocyte")
p<-DotPlot(BM.combined.sct, cluster.idents = F,features = c("CD34","ENG","ITGAM","CD3G","CD19","FCGR3A","CEACAM8","ITGA2B","GP1BA"), assay = "SCT",scale.min=0,scale.max=60)+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+ scale_colour_gradient2(low = "white", high = "#107dac")+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p1<-DotPlot(BM.combined.sct,cluster.idents = F, features = c("CD34.1","CD105","CD11b","CD3","CD19.1","CD16","CD66b","CD41",'CD42b'), assay = "ADT2",scale.min=0,scale.max=100)+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+ scale_colour_gradient2(low = "white", high = "#BF3F3F")+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
multiplot(p,p1, cols = 2)

Idents(BM.combined.sct) = "manualAnno_specific"
levels(BM.combined.sct) <- c("Differentiating Stem Cell","CD34","Stromal","EPC","Monocyte","DC","NK","gdT","CD4","CD8","MAIT","B cell","Plasmablast","Neutrophil","Megakaryocyte")
pp<-DotPlot(BM.combined.sct, cluster.idents = F,features = c(c("CD34","PDGFRA","ENG","TFRC","CD14","ITGAM","ITGAX","CD1C","FCGR3A","NCAM1","TRGV9","CD4","CD8A","KLRB1","CD19","MS4A1","CEACAM8","ITGA2B","GP1BA")), assay = "SCT",scale.min=0,scale.max=60)+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+ scale_colour_gradient2(low = "white", high = "#107dac")+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
pp1<-DotPlot(BM.combined.sct,cluster.idents = F, features = c("CD34.1","CD140a-PDGFRa","CD105","CD71","CD14.1","CD11b","CD11c","CD1c","CD16","CD56","TCR-Vy9","CD4.1","CD8","CD161","CD19.1","CD20","CD66b","CD41"), assay = "ADT2",scale.min=0,scale.max=100)+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+ scale_colour_gradient2(low = "white", high = "#BF3F3F",midpoint = 0)+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
multiplot(pp,pp1, cols = 2)

Idents(BM.combined.sct) = "manualAnno_specific2"
levels(BM.combined.sct) <- rev(c("Stromal","EPC","MEP","HSC","LMPP","GMP","Neutrophil","Classical Mono","Intermediate Mono","Non-classical Mono","mDC","pDC","CD56dim NK","CD56bright NK","Vd1 gdT","Vd2 gdT","CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL","CD4 Treg","CD8 Naive","CD8 TCM","CD8 TEM","CD8 TEMRA","MAIT","Naive B","Memory B","Plasmablast","Megakaryocyte","Differentiating Stem Cell"))
ppp<-DotPlot(BM.combined.sct, cluster.idents = F,features = c("CXCL12","ENG","TFRC","HPGDS","GATA1","GATA2","CD34","KIT","CRHBP","AVP","NPR3","CD38","NPW","PRTN3","ELANE","AZU1","CEACAM8","CD14","ITGAM","ITGAX","CD1C","CLEC4C","FCGR3A","NCAM1","TRGV9","TRDV1","TRDV2","FOXP3","CTLA4","CD4","CD8A","CCR7","SELL","CD28","KLRG1","KLRB1","CD19","MS4A1","IGHM","IGHD","IGHG2","CD27","ITGA2B","GP1BA"), assay = "SCT",scale.min=0,scale.max=60)+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+ scale_colour_gradient2(low = "white", high = "#107dac")+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ppp1<-DotPlot(BM.combined.sct,cluster.idents = F, features = c("CD105","CD71","CD34.1","CD90-Thy1","CD49f","CD38.1","CD66b","HLA-DR","CD14.1","CD11b","CD11c","CD1c","CD303-BDCA2","CD16","CD56","TCR-Vy9","TCR-Vd2","CD25","CD4.1","CD8","CD45RA","CD45RO","CD62L","CD28.1","KLRG1-MAFA","CD161","CD19.1","CD20","IgD","IgM","CD27.1","CD41","CD42b"), assay = "ADT2",scale.min=0,scale.max=100)+geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+ scale_colour_gradient2(low = "white", high = "#BF3F3F")+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
multiplot(ppp,ppp1, cols = 2)


SaveSeuratRds(BM.combined.sct, "BM.combined.sct_2025.rds")


####################################              IMMUNE FREQUENCY ANALYSIS               ####################################
BM.combined.sct <- LoadSeuratRds("BM.combined.sct_2025.rds")
BM.combined.sct$SampleInfo_Type_Time <- paste(BM.combined.sct$Patient,BM.combined.sct$SampleType,BM.combined.sct$Time,sep="-")
BM.combined.sct$SampleInfo_Type_Time<- factor(BM.combined.sct$SampleInfo_Type_Time, levels = c("24780-BM-Scr","24780-BM-PC2","122759-BM-Scr","122759-BM-PC2","113446-BM-Scr","113446-BM-PC2","122739-BM-Scr","122739-BM-PC2","121535-BM-Scr","121535-BM-PC2","121474-BM-Scr","121474-BM-PC2","119837-BM-Scr","119837-BM-PC2","121015-BM-Scr","121015-BM-PC2","122976-BM-Scr","122976-BM-PC2","114152-BM-Scr","114152-BM-PC2","119325-BM-Scr","119325-BM-PC2","118877-BM-Scr","118877-BM-PC2","115585-BM-Scr","115585-BM-PC2","24780-PBMC-Scr","24780-PBMC-C2D1","24780-PBMC-C2D8","121474-PBMC-Scr","121474-PBMC-C1D8","121474-PBMC-C2D1","121015-PBMC-Scr","121015-PBMC-C1D8","121015-PBMC-C2D1","118877-PBMC-Scr","118877-PBMC-C1D8","118877-PBMC-C2D1"))

BM.combined.sct@meta.data$PatientTime <- paste(BM.combined.sct@meta.data$Patient,BM.combined.sct@meta.data$Time2, sep="_")
BM.combined.sct@meta.data$TimeSurvival <- paste(BM.combined.sct@meta.data$Survival,BM.combined.sct@meta.data$Time2,sep = "_")
BM.combined.sct@meta.data$Time2 <- factor(BM.combined.sct@meta.data$Time2, levels = c("Pre","C1D8","C2D8","Post"))

BM.combined.sct.BM <- subset(BM.combined.sct, subset = SampleType == "BM")
BM.combined.sct.PBMC <- subset(BM.combined.sct, subset = SampleType == "PBMC")



t = data.frame()
for (i in as.character(unique(BM.combined.sct$SampleInfo_Type_Time))) 
{d<-data.frame(table(BM.combined.sct@meta.data[BM.combined.sct$SampleInfo_Type_Time == i,]$manualAnno_specific))
for (cell in unique(BM.combined.sct$manualAnno_specific)) {if (!(cell %in% d$Var1)==T) {d<-rbind(d,data.frame(Var1=cell,Freq=0))}}
d$TotalCells <- nrow(BM.combined.sct@meta.data[BM.combined.sct$SampleInfo_Type_Time == i,])
d$SampleInfo_Type_Time <- i
d$Patient <- unique(BM.combined.sct@meta.data[BM.combined.sct@meta.data$SampleInfo_Type_Time %in% i,]$Patient)
d$Time2<- unique(BM.combined.sct@meta.data[BM.combined.sct@meta.data$SampleInfo_Type_Time %in% i,]$Time2)
d$SampleType<- unique(BM.combined.sct@meta.data[BM.combined.sct@meta.data$SampleInfo_Type_Time %in% i,]$SampleType)
d$Survival<- unique(BM.combined.sct@meta.data[BM.combined.sct@meta.data$SampleInfo_Type_Time %in% i,]$Survival)
d$Gender<- unique(BM.combined.sct@meta.data[BM.combined.sct@meta.data$SampleInfo_Type_Time %in% i,]$Gender)
d$Race<- unique(BM.combined.sct@meta.data[BM.combined.sct@meta.data$SampleInfo_Type_Time %in% i,]$Race)
d$type <- "AnnoSpecific"
d$Time<- unique(BM.combined.sct@meta.data[BM.combined.sct@meta.data$SampleInfo_Type_Time %in% i,]$Time)

t<-rbind(t,d)
}

t$SampleInfo_Type_Time <-factor(t$SampleInfo_Type_Time,levels = c("24780-BM-Scr","24780-BM-PC2","122759-BM-Scr","122759-BM-PC2","113446-BM-Scr","113446-BM-PC2","122739-BM-Scr","122739-BM-PC2","121535-BM-Scr","121535-BM-PC2","121474-BM-Scr","121474-BM-PC2","119837-BM-Scr","119837-BM-PC2","121015-BM-Scr","121015-BM-PC2","122976-BM-Scr","122976-BM-PC2","114152-BM-Scr","114152-BM-PC2","119325-BM-Scr","119325-BM-PC2","118877-BM-Scr","118877-BM-PC2","115585-BM-Scr","115585-BM-PC2","24780-PBMC-Scr","24780-PBMC-C2D1","24780-PBMC-C2D8","121474-PBMC-Scr","121474-PBMC-C1D8","121474-PBMC-C2D1","121015-PBMC-Scr","121015-PBMC-C1D8","121015-PBMC-C2D1","118877-PBMC-Scr","118877-PBMC-C1D8","118877-PBMC-C2D1"))
t$Var1 <- factor(t$Var1,levels = c("CD34","EPC","Monocyte","DC","B cell","Plasmablast","gdT","CD4","CD8","MAIT","NK","Stromal","Megakaryocyte","Neutrophil","Differentiating Stem Cell"))
t <- t[order(t$SampleInfo_Type_Time),]
t$percentage <- round(t$Freq/t$TotalCells,3)

totalcell_count <- unique(t$TotalCells)
median(totalcell_count) #2242.5
mean(totalcell_count) #2251
min(totalcell_count) #78
max(totalcell_count) # 5329

t_BMtest <- t[t$SampleType == "BM",]
totalcell_count_BM <- unique(t_BMtest$TotalCells)
median(totalcell_count_BM) #2242.5
mean(totalcell_count_BM) #2251
min(totalcell_count_BM) #78
max(totalcell_count_BM) # 5329

t_PBMCtest <- t[t$SampleType == "PBMC",]
totalcell_count_PBMC <- unique(t_PBMCtest$TotalCells) #24051
median(totalcell_count_PBMC) #2065
mean(totalcell_count_PBMC) #2004
min(totalcell_count_PBMC) #874
max(totalcell_count_PBMC) # 4336

pp<-ggplot(t, aes(x=SampleInfo_Type_Time, y=Freq, fill=Var1)) +geom_bar(stat="identity",colour = "black", position="fill")+ggtitle("Immune Composition")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Patient Samples", y = "Frequency")+scale_fill_manual(values = c('#ff0000', '#dbac98', '#ffe700', '#c27ba0', '#6565bf','#660066', '#e7b416', '#99cc33','#02a9f7','#7df9ff','#ff7400','#7FFFCC', 'lightgrey','darkgrey','#ffdfba'))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
pp
p<-ggplot(t, aes(x=SampleInfo_Type_Time, y=Freq, fill=Var1)) +geom_bar(stat="identity",colour = "black", position="stack")+ggtitle("Immune Composition")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Patient Samples", y = "Frequency")+scale_fill_manual(values = c('#ff0000', '#dbac98', '#ffe700', '#c27ba0', '#6565bf','#660066', '#e7b416', '#99cc33','#02a9f7','#7df9ff','#ff7400','#7FFFCC', 'lightgrey','darkgrey','#ffdfba'))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
p

patchwork::wrap_plots(p,pp, ncol = 1)


## Plot heatmap of % change of immune frequency by patient
t_BM <- t[t$SampleType == "BM",c("Patient","Time2","Var1","Freq","TotalCells","percentage")]
t_BM_diff = data.frame()
for (i in unique(t_BM$Patient))
{d <- t_BM[t_BM$Patient == i & t_BM$Time2 == "Pre", c("Patient","Var1","percentage")]
d2 <- t_BM[t_BM$Patient == i & t_BM$Time2 == "Post",c("Patient","Var1","percentage")] 
d3 <-merge(d,d2, by = c("Patient","Var1"))
t_BM_diff<-rbind(t_BM_diff,d3)
}

t_BM_diff$DiffPerc <- t_BM_diff$percentage.y-t_BM_diff$percentage.x
t_BM_diff$Patient <- factor(t_BM_diff$Patient, levels = c("24780","122759","113446","122739","121535","121474","119837","121015","122976","114152","119325","118877","115585"))

ggplot(t_BM_diff, aes(y = Var1, x = Patient, fill = DiffPerc))+geom_tile()+scale_y_discrete(limits=rev)+ scale_fill_gradient2(high = "#bf9000", low = "#3E8245",limits = c(-0.5,0.5)) +theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))


t_BM_diff[t_BM_diff$DiffPerc >0.20,]$DiffPerc <- 0.20 #for visualization purposes, change anything > 0.2 to 0.2
t_BM_diff[t_BM_diff$DiffPerc < -0.20,]$DiffPerc <- -0.20 #for visualization purposes, change anything < -0.2 to -0.2

ggplot(t_BM_diff, aes(y = Var1, x = Patient, fill = DiffPerc))+geom_tile()+scale_y_discrete(limits=rev)+ scale_fill_gradient2(high = "#bf9000", low = "#3E8245",limits = c(-0.2,0.2)) +theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))




## Plot heatmap of % change of immune frequency by patient (PBMC)
t_PBMC <- t[t$SampleType == "PBMC",c("Patient","Time","Var1","Freq","TotalCells","percentage")]
t_PBMC_diff = data.frame()
for (i in 24780)
{d <- t_PBMC[t_PBMC$Patient == i & t_PBMC$Time == "Scr", c("Patient","Var1","percentage")]
d3 <- t_PBMC[t_PBMC$Patient == i & t_PBMC$Time == "C2D1",c("Patient","Var1","percentage")] 
d4 <- t_PBMC[t_PBMC$Patient == i & t_PBMC$Time == "C2D8",c("Patient","Var1","percentage")] 
df_list <- list(d,d3, d4)
d5 <- Reduce(function(x, y) merge(x, y,  by = c("Patient","Var1")), df_list, accumulate=FALSE)
t_PBMC_diff<-rbind(t_PBMC_diff,d5)
}
for (i in c(121474,121015,118877))
{d <- t_PBMC[t_PBMC$Patient == i & t_PBMC$Time == "Scr", c("Patient","Var1","percentage")]
d2 <- t_PBMC[t_PBMC$Patient == i & t_PBMC$Time == "C1D8",c("Patient","Var1","percentage")] 
d3 <- t_PBMC[t_PBMC$Patient == i & t_PBMC$Time == "C2D1",c("Patient","Var1","percentage")] 
df_list <- list(d,d2, d3)
d5 <- Reduce(function(x, y) merge(x, y,  by = c("Patient","Var1")), df_list, accumulate=FALSE)
t_PBMC_diff<-rbind(t_PBMC_diff,d5)
}

t_PBMC_diff$DiffPerc <- t_PBMC_diff$percentage.y-t_PBMC_diff$percentage.x
t_PBMC_diff$DiffPerc1 <- t_PBMC_diff$percentage-t_PBMC_diff$percentage.x
t_PBMC_diff$DiffPerc2 <- t_PBMC_diff$percentage-t_PBMC_diff$percentage.y

t_PBMC_diff$Patient <- factor(t_PBMC_diff$Patient, levels = c("24780","122759","113446","122739","121535","121474","119837","121015","122976","114152","119325","118877","115585"))

a<-ggplot(t_PBMC_diff, aes(y = Var1, x = Patient, fill = DiffPerc))+geom_tile()+scale_y_discrete(limits=rev)+ scale_fill_gradient2(high = "#bf9000", low = "#3E8245",limits = c(-0.5,0.5)) +theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
b<-ggplot(t_PBMC_diff, aes(y = Var1, x = Patient, fill = DiffPerc1))+geom_tile()+scale_y_discrete(limits=rev)+ scale_fill_gradient2(high = "#bf9000", low = "#3E8245",limits = c(-0.5,0.5)) +theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
c<-ggplot(t_PBMC_diff, aes(y = Var1, x = Patient, fill = DiffPerc2))+geom_tile()+scale_y_discrete(limits=rev)+ scale_fill_gradient2(high = "#bf9000", low = "#3E8245",limits = c(-0.5,0.5)) +theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))


t_PBMC_diff[t_PBMC_diff$DiffPerc >0.15,]$DiffPerc <- 0.15 #for visualization purposes, change anything > 0.15 to 0.15
t_PBMC_diff[t_PBMC_diff$DiffPerc < -0.15,]$DiffPerc <- -0.15 #for visualization purposes, change anything < -0.15 to -0.15
t_PBMC_diff[t_PBMC_diff$DiffPerc1 >0.15,]$DiffPerc1 <- 0.15 #for visualization purposes, change anything > 0.15 to 0.15
t_PBMC_diff[t_PBMC_diff$DiffPerc1 < -0.15,]$DiffPerc1 <- -0.15 #for visualization purposes, change anything < -0.15 to -0.15
t_PBMC_diff[t_PBMC_diff$DiffPerc2 >0.15,]$DiffPerc2 <- 0.15 #for visualization purposes, change anything > 0.15 to 0.15
t_PBMC_diff[t_PBMC_diff$DiffPerc2 < -0.15,]$DiffPerc2 <- -0.15 #for visualization purposes, change anything < -0.15 to -0.15


a<-ggplot(t_PBMC_diff, aes(y = Var1, x = Patient, fill = DiffPerc))+geom_tile()+scale_y_discrete(limits=rev)+ scale_fill_gradient2(high = "#bf9000", low = "#3E8245",limits = c(-0.15,0.15)) +theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
b<-ggplot(t_PBMC_diff, aes(y = Var1, x = Patient, fill = DiffPerc1))+geom_tile()+scale_y_discrete(limits=rev)+ scale_fill_gradient2(high = "#bf9000", low = "#3E8245",limits = c(-0.15,0.15)) +theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
c<-ggplot(t_PBMC_diff, aes(y = Var1, x = Patient, fill = DiffPerc2))+geom_tile()+scale_y_discrete(limits=rev)+ scale_fill_gradient2(high = "#bf9000", low = "#3E8245",limits = c(-0.15,0.15)) +theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))

multiplot(a,b,c)




#### Specific Anno Immune % calculation
t2 = data.frame()
for (i in unique(BM.combined.sct$SampleInfo_Type_Time)) 
{d<-data.frame(table(BM.combined.sct@meta.data[BM.combined.sct$SampleInfo_Type_Time == i,]$manualAnno_specific2))
for (cell in unique(BM.combined.sct$manualAnno_specific2)) {if (!(cell %in% d$Var1)==T) {d<-rbind(d,data.frame(Var1=cell,Freq=0))}}
d$TotalCells <- nrow(BM.combined.sct@meta.data[BM.combined.sct$SampleInfo_Type_Time == i,])
d$SampleInfo_Type_Time <- i
d$Patient <- unique(BM.combined.sct@meta.data[BM.combined.sct$SampleInfo_Type_Time == i,]$Patient)
d$Time2<- unique(BM.combined.sct@meta.data[BM.combined.sct$SampleInfo_Type_Time == i,]$Time2)
d$SampleType<- unique(BM.combined.sct@meta.data[BM.combined.sct$SampleInfo_Type_Time == i,]$SampleType)
d$Survival<- unique(BM.combined.sct@meta.data[BM.combined.sct$SampleInfo_Type_Time == i,]$Survival)
d$Gender<- unique(BM.combined.sct@meta.data[BM.combined.sct$SampleInfo_Type_Time == i,]$Gender)
d$Race<- unique(BM.combined.sct@meta.data[BM.combined.sct$SampleInfo_Type_Time == i,]$Race)
d$type <- "AnnoSpecific"
d$Time<- unique(BM.combined.sct@meta.data[BM.combined.sct$SampleInfo_Type_Time == i,]$Time)
t2<-rbind(t2,d)
}

t2$SampleInfo_Type_Time <-factor(t2$SampleInfo_Type_Time,levels = c("24780-BM-Scr","24780-BM-PC2","122759-BM-Scr","122759-BM-PC2","113446-BM-Scr","113446-BM-PC2","122739-BM-Scr","122739-BM-PC2","121535-BM-Scr","121535-BM-PC2","121474-BM-Scr","121474-BM-PC2","119837-BM-Scr","119837-BM-PC2","121015-BM-Scr","121015-BM-PC2","122976-BM-Scr","122976-BM-PC2","114152-BM-Scr","114152-BM-PC2","119325-BM-Scr","119325-BM-PC2","118877-BM-Scr","118877-BM-PC2","115585-BM-Scr","115585-BM-PC2","24780-PBMC-Scr","24780-PBMC-C2D1","24780-PBMC-C2D8","121474-PBMC-Scr","121474-PBMC-C1D8","121474-PBMC-C2D1","121015-PBMC-Scr","121015-PBMC-C1D8","121015-PBMC-C2D1","118877-PBMC-Scr","118877-PBMC-C1D8","118877-PBMC-C2D1"))
t2$Var1 <- factor(t2$Var1,levels = c("EPC","MEP","HSC","LMPP","GMP","Classical Mono","Intermediate Mono","Non-classical Mono","pDC","mDC","Naive B","Memory B","Plasmablast","Vd1 gdT","Vd2 gdT","CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL","CD4 Treg","CD8 Naive","CD8 TCM","CD8 TEM","CD8 TEMRA","MAIT","CD56dim NK","CD56bright NK",'Stromal', "Megakaryocyte","Neutrophil","Differentiating Stem Cell"))
t2 <- t2[order(t2$SampleInfo_Type_Time),]
t2$percentage <- round(t2$Freq/t2$TotalCells,3)

pp<-ggplot(t2, aes(x=SampleInfo_Type_Time, y=Freq, fill=Var1)) +geom_bar(stat="identity",colour = "black", position="fill")+ggtitle("Immune Composition")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Patient Samples", y = "Frequency")+scale_fill_manual(values = c('#dbac98','#a70000','#ff0000','#bca9bd','#ffbaf2','#d0b783','yellow','#876127','#326ada','#c27ba0','#efbbff','#6565bf','#660066','#ffd600','#e7b416',"#339900","#99cc33","#0c7741","#8cff32","#e5c3c6",'#bae1ff','#b3cde0','#005b96','#5bc0de','#7df9ff','#ff9966','#ff4f00','#7FFFCC','lightgrey','darkgrey','#ffdfba'))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
pp
p<-ggplot(t2, aes(x=SampleInfo_Type_Time, y=Freq, fill=Var1)) +geom_bar(stat="identity",colour = "black", position="stack")+ggtitle("Immune Composition")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Patient Samples", y = "Frequency")+scale_fill_manual(values = c('#dbac98','#a70000','#ff0000','#bca9bd','#ffbaf2','#d0b783','yellow','#876127','#326ada','#c27ba0','#efbbff','#6565bf','#660066','#ffd600','#e7b416',"#339900","#99cc33","#0c7741","#8cff32","#e5c3c6",'#bae1ff','#b3cde0','#005b96','#5bc0de','#7df9ff','#ff9966','#ff4f00','#7FFFCC','lightgrey','darkgrey','#ffdfba'))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
p

patchwork::wrap_plots(p,pp, ncol = 1)


## Plot heatmap of % change of immune frequency by patient
t2_BM <- t2[t2$SampleType == "BM",c("Patient","Time2","Var1","Freq","TotalCells","percentage")]
t2_BM_diff = data.frame()
for (i in unique(t2_BM$Patient))
{d <- t2_BM[t2_BM$Patient == i & t2_BM$Time2 == "Pre", c("Patient","Var1","percentage")]
d2 <- t2_BM[t2_BM$Patient == i & t2_BM$Time2 == "Post",c("Patient","Var1","percentage")] 
d3 <-merge(d,d2, by = c("Patient","Var1"))
t2_BM_diff<-rbind(t2_BM_diff,d3)
}

t2_BM_diff$DiffPerc <- t2_BM_diff$percentage.y-t2_BM_diff$percentage.x
t2_BM_diff$Patient <- factor(t2_BM_diff$Patient, levels = c("24780","122759","113446","122739","121535","121474","119837","121015","122976","114152","119325","118877","115585"))

ggplot(t2_BM_diff, aes(y = Var1, x = Patient, fill = DiffPerc))+geom_tile()+scale_y_discrete(limits=rev)+ scale_fill_gradient2(high = "#bf9000", low = "#3E8245",limits = c(-0.5,0.5)) +theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))

t2_BM_diff[t2_BM_diff$DiffPerc >0.10,]$DiffPerc <- 0.10 #for visualization purposes, change anything > 0.2 to 0.2
t2_BM_diff[t2_BM_diff$DiffPerc < -0.10,]$DiffPerc <- -0.10 #for visualization purposes, change anything < -0.2 to -0.2

ggplot(t2_BM_diff, aes(y = Var1, x = Patient, fill = DiffPerc))+geom_tile()+scale_y_discrete(limits=rev)+ scale_fill_gradient2(high = "#bf9000", low = "#3E8245",limits = c(-0.10,0.10)) +theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))


t2_BM_diff$percentage.x1 <- t2_BM_diff$percentage.x
t2_BM_diff[t2_BM_diff$percentage.x1 >0.30,]$percentage.x1 <- 0.30 #for visualization purposes, change anything > 0.1 to 0.1
t2_BM_diff$percentage.y1 <- t2_BM_diff$percentage.y
t2_BM_diff[t2_BM_diff$percentage.y1 >0.30,]$percentage.y1 <- 0.30 #for visualization purposes, change anything > 0.1 to 0.1


aa<-ggplot(t2_BM_diff, aes(y = Var1, x = Patient, fill = percentage.x1))+geom_tile()+scale_y_discrete(limits=rev)+ scale_fill_gradient2(high = "#777696", low = "white",limits = c(0,0.3)) +theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
bb<-ggplot(t2_BM_diff, aes(y = Var1, x = Patient, fill = percentage.y1))+geom_tile()+scale_y_discrete(limits=rev)+ scale_fill_gradient2(high = "#777696", low = "white",limits = c(0,0.3)) +theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
multiplot(aa,bb, cols = 2)





## Plot heatmap of % change of immune frequency by patient (PBMC)
t2_PBMC <- t2[t2$SampleType == "PBMC",c("Patient","Time","Var1","Freq","TotalCells","percentage")]
t2_PBMC_diff = data.frame()
for (i in 24780)
{d <- t2_PBMC[t2_PBMC$Patient == i & t2_PBMC$Time == "Scr", c("Patient","Var1","percentage")]
d3 <- t2_PBMC[t2_PBMC$Patient == i & t2_PBMC$Time == "C2D1",c("Patient","Var1","percentage")] 
d4 <- t2_PBMC[t2_PBMC$Patient == i & t2_PBMC$Time == "C2D8",c("Patient","Var1","percentage")] 
df_list <- list(d,d3, d4)
d5 <- Reduce(function(x, y) merge(x, y,  by = c("Patient","Var1")), df_list, accumulate=FALSE)
t2_PBMC_diff<-rbind(t2_PBMC_diff,d5)
}
for (i in c(121474,121015,118877))
{d <- t2_PBMC[t2_PBMC$Patient == i & t2_PBMC$Time == "Scr", c("Patient","Var1","percentage")]
d2 <- t2_PBMC[t2_PBMC$Patient == i & t2_PBMC$Time == "C1D8",c("Patient","Var1","percentage")] 
d3 <- t2_PBMC[t2_PBMC$Patient == i & t2_PBMC$Time == "C2D1",c("Patient","Var1","percentage")] 
df_list <- list(d,d2, d3)
d5 <- Reduce(function(x, y) merge(x, y,  by = c("Patient","Var1")), df_list, accumulate=FALSE)
t2_PBMC_diff<-rbind(t2_PBMC_diff,d5)
}

t2_PBMC_diff$DiffPerc <- t2_PBMC_diff$percentage.y-t2_PBMC_diff$percentage.x
t2_PBMC_diff$DiffPerc1 <- t2_PBMC_diff$percentage-t2_PBMC_diff$percentage.x
t2_PBMC_diff$DiffPerc2 <- t2_PBMC_diff$percentage-t2_PBMC_diff$percentage.y

t2_PBMC_diff$Patient <- factor(t2_PBMC_diff$Patient, levels = c("24780","122759","113446","122739","121535","121474","119837","121015","122976","114152","119325","118877","115585"))

a<-ggplot(t2_PBMC_diff, aes(y = Var1, x = Patient, fill = DiffPerc))+geom_tile()+scale_y_discrete(limits=rev)+ scale_fill_gradient2(high = "#bf9000", low = "#3E8245",limits = c(-0.5,0.5)) +theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
b<-ggplot(t2_PBMC_diff, aes(y = Var1, x = Patient, fill = DiffPerc1))+geom_tile()+scale_y_discrete(limits=rev)+ scale_fill_gradient2(high = "#bf9000", low = "#3E8245",limits = c(-0.5,0.5)) +theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
c<-ggplot(t2_PBMC_diff, aes(y = Var1, x = Patient, fill = DiffPerc2))+geom_tile()+scale_y_discrete(limits=rev)+ scale_fill_gradient2(high = "#bf9000", low = "#3E8245",limits = c(-0.5,0.5)) +theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))


t2_PBMC_diff[t2_PBMC_diff$DiffPerc >0.10,]$DiffPerc <- 0.10 #for visualization purposes, change anything > 0.1 to 0.1
t2_PBMC_diff[t2_PBMC_diff$DiffPerc < -0.10,]$DiffPerc <- -0.10 #for visualization purposes, change anything < -0.1 to -0.1
t2_PBMC_diff[t2_PBMC_diff$DiffPerc1 >0.10,]$DiffPerc1 <- 0.10 #for visualization purposes, change anything > 0.1 to 0.1
t2_PBMC_diff[t2_PBMC_diff$DiffPerc1 < -0.10,]$DiffPerc1 <- -0.10 #for visualization purposes, change anything < -0.1 to -0.1
t2_PBMC_diff[t2_PBMC_diff$DiffPerc2 >0.10,]$DiffPerc2 <- 0.10 #for visualization purposes, change anything > 0.1 to 0.1
t2_PBMC_diff[t2_PBMC_diff$DiffPerc2 < -0.10,]$DiffPerc2 <- -0.10 #for visualization purposes, change anything < -0.1 to -0.1


a<-ggplot(t2_PBMC_diff, aes(y = Var1, x = Patient, fill = DiffPerc))+geom_tile()+scale_y_discrete(limits=rev)+ scale_fill_gradient2(high = "#bf9000", low = "#3E8245",limits = c(-0.1,0.1)) +theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
b<-ggplot(t2_PBMC_diff, aes(y = Var1, x = Patient, fill = DiffPerc1))+geom_tile()+scale_y_discrete(limits=rev)+ scale_fill_gradient2(high = "#bf9000", low = "#3E8245",limits = c(-0.1,0.1)) +theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
c<-ggplot(t2_PBMC_diff, aes(y = Var1, x = Patient, fill = DiffPerc2))+geom_tile()+scale_y_discrete(limits=rev)+ scale_fill_gradient2(high = "#bf9000", low = "#3E8245",limits = c(-0.1,0.1)) +theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))

multiplot(a,b,c, cols = 3)

t2_PBMC_diff$percentage.x1 <- t2_PBMC_diff$percentage.x
t2_PBMC_diff[t2_PBMC_diff$percentage.x1 >0.30,]$percentage.x1 <- 0.30 #for visualization purposes, change anything > 0.1 to 0.1
t2_PBMC_diff$percentage.y1 <- t2_PBMC_diff$percentage.y
t2_PBMC_diff[t2_PBMC_diff$percentage.y1 >0.30,]$percentage.y1 <- 0.30 #for visualization purposes, change anything > 0.1 to 0.1
t2_PBMC_diff$percentage.1 <- t2_PBMC_diff$percentage
t2_PBMC_diff[t2_PBMC_diff$percentage.1 >0.30,]$percentage.1 <- 0.30 #for visualization purposes, change anything > 0.1 to 0.1


aa<-ggplot(t2_PBMC_diff, aes(y = Var1, x = Patient, fill = percentage.x1))+geom_tile()+scale_y_discrete(limits=rev)+ scale_fill_gradient2(high = "#740001", low = "white",limits = c(0,0.3)) +theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
bb<-ggplot(t2_PBMC_diff, aes(y = Var1, x = Patient, fill = percentage.y1))+geom_tile()+scale_y_discrete(limits=rev)+ scale_fill_gradient2(high = "#740001", low = "white",limits = c(0,0.3)) +theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
cc<-ggplot(t2_PBMC_diff, aes(y = Var1, x = Patient, fill = percentage.1))+geom_tile()+scale_y_discrete(limits=rev)+ scale_fill_gradient2(high = "#740001", low = "white",limits = c(0,0.3)) +theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
multiplot(aa,bb,cc, cols = 3)


############### COMPARE CITE-SEQ ANNOTATION WITH ORTHOGONAL BM FACS EXPERIMENT #############
t_FACS <- t[t$Patient %in% c(119325,119837,121015,121535,122976) & t$Time2 == "Pre" & t$SampleType == "BM",]
#specific annotation
t_FACS_specific <- t_FACS
t_FACS_specific$Var1 <- as.character(t_FACS_specific$Var1)
t_FACS_specific[t_FACS_specific$Var1 == "Monocyte",]$Var1 <- "Myeloid"
t_FACS_specific[t_FACS_specific$Var1 == "DC",]$Var1 <- "Myeloid"
t_FACS_specific[t_FACS_specific$Var1 == "Plasmablast",]$Var1 <- "B cell"
t_FACS_specific[t_FACS_specific$Var1 == "MAIT",]$Var1 <- "Others"
t_FACS_specific[t_FACS_specific$Var1 == "EPC",]$Var1 <- "Others"
t_FACS_specific[t_FACS_specific$Var1 == "Differentiating Stem Cell",]$Var1 <- "Others"
t_FACS_specific[t_FACS_specific$Var1 == "Stromal",]$Var1 <- "Others"
t_FACS_specific[t_FACS_specific$Var1 == "Neutrophil",]$Var1 <- "Others"
t_FACS_specific[t_FACS_specific$Var1 == "Megakaryocyte",]$Var1 <- "Others"

t_FACS_specific_agg <- aggregate(Freq ~ Var1 + Patient + TotalCells, t_FACS_specific[,c(1,2,3,5)], FUN = sum)
t_FACS_specific_agg$percentage <- round(t_FACS_specific_agg$Freq/t_FACS_specific_agg$TotalCells,2)
t_FACS_specific_agg$Var1 <- factor(t_FACS_specific_agg$Var1,levels = c("CD34","Myeloid", "B cell","gdT", "CD4", "CD8","NK","Others"))
t_FACS_specific_agg$Patient <- factor(t_FACS_specific_agg$Patient,levels = c("122976","121535", "121015","119837", "119325"))


p1<-ggplot(t_FACS_specific_agg, aes(x=Patient, y=percentage, fill=Var1)) +geom_bar(stat="identity",colour = "black", position="stack")+ggtitle("CITE-seq BM Immune Composition")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Patient Samples", y = "Frequency")+scale_fill_manual(values = c('#ff0000', '#ffe700', '#6565bf', '#e7b416',"#99cc33",'#02a9f7','orange','lightgrey'))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
p1

FACS_specific <- read.table("2026-01-09_MDS_NewGating_Counts.txt", sep = "\t", header = T, stringsAsFactors = F)
FACS_specific <- FACS_specific[c(1:5),c(1,4:14)]
colnames(FACS_specific) <- c("Patient", "TotalCells","CD34","EPC","Myeloid", "NK", "CD11bnegCD56neg","gdT","CD4","CD8","NKT-like", "B cell")
FACS_specific$Others <- FACS_specific$TotalCells-(rowSums(FACS_specific[,c(3,5,6,8,9,10,12)])) #include CD11bnegCD56neg and EPC and NKT-like in Others
FACS_specific <- FACS_specific[,c(1:3,5,6,8,9,10,12,13)]
FACS_specific<-reshape2::melt(FACS_specific, id = c("Patient","TotalCells"))
FACS_specific$percentage <- round(FACS_specific$value/FACS_specific$TotalCells,2)
FACS_specific$variable <- factor(FACS_specific$variable,levels =  c("CD34","Myeloid", "B cell","gdT", "CD4", "CD8","NK","Others"))
FACS_specific$Patient <- factor(FACS_specific$Patient,levels = c("122976.fcs","121535.fcs", "121015.fcs","119837.fcs", "119325.fcs"))

pp1<-ggplot(FACS_specific, aes(x=Patient, y=percentage, fill=variable)) +geom_bar(stat="identity",colour = "black", position="stack")+ggtitle("CITE-seq BM Immune Composition")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Patient Samples", y = "Frequency")+scale_fill_manual(values = c('#ff0000', '#ffe700', '#6565bf', '#e7b416',"#99cc33",'#02a9f7','orange','lightgrey'))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
pp1

multiplot(p1,pp1,cols = 2)             


##### PERFORM WELCH T-TEST FOR DIFFERENTIAL ABUNDANCE ##################
#Short_Pre vs Short_Post
for (i in unique(t$Var1)) {
  print(i)
  print(t.test(t[t$Var1 == i & t$Time2 == "Pre" & t$Survival == "Short" & t$SampleType == "BM",]$percentage, t[t$Var1 == i & t$Time2 == "Post" & t$Survival == "Short" & t$SampleType == "BM",]$percentage))}


#Long_Pre vs Long_Post
for (i in unique(t$Var1)) {
  print(i)
  print(t.test(t[t$Var1 == i & t$Time2 == "Pre" & t$Survival == "Long" & t$SampleType == "BM",]$percentage, t[t$Var1 == i & t$Time2 == "Post" & t$Survival == "Long" & t$SampleType == "BM",]$percentage))}

[1] "CD34"

Welch Two Sample t-test

data:  t[t$Var1 == i & t$Time2 == "Pre" & t$Survival == "Long" & t$SampleType == "BM", ]$percentage and t[t$Var1 == i & t$Time2 == "Post" & t$Survival == "Long" & t$SampleType == "BM", ]$percentage
t = 2.1862, df = 12.576, p-value = 0.04836
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  0.001050278 0.248149722
sample estimates:
  mean of x mean of y 
0.2290    0.1044 




#Short_Pre vs Long_Pre
for (i in unique(t$Var1)) {
  print(i)
  print(t.test(t[t$Var1 == i & t$Time2 == "Pre" & t$Survival == "Short" & t$SampleType == "BM",]$percentage, t[t$Var1 == i & t$Time2 == "Pre" & t$Survival == "Long" & t$SampleType == "BM",]$percentage))}

[1] "Monocyte"

Welch Two Sample t-test

data:  t[t$Var1 == i & t$Time2 == "Pre" & t$Survival == "Short" & t$SampleType == "BM", ]$percentage and t[t$Var1 == i & t$Time2 == "Pre" & t$Survival == "Long" & t$SampleType == "BM", ]$percentage
t = -2.2455, df = 9.717, p-value = 0.04929
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -0.2053436331 -0.0003897003
sample estimates:
  mean of x  mean of y 
0.01633333 0.11920000 


[1] "CD4"

Welch Two Sample t-test

data:  t[t$Var1 == i & t$Time2 == "Pre" & t$Survival == "Short" & t$SampleType == "BM", ]$percentage and t[t$Var1 == i & t$Time2 == "Pre" & t$Survival == "Long" & t$SampleType == "BM", ]$percentage
t = 2.9654, df = 3.5106, p-value = 0.04868
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  0.00307182 0.61572818
sample estimates:
  mean of x mean of y 
0.5360    0.2266 



#Short_Post vs Long_Post
for (i in unique(t$Var1)) {
  print(i)
  print(t.test(t[t$Var1 == i & t$Time2 == "Post" & t$Survival == "Short" & t$SampleType == "BM",]$percentage, t[t$Var1 == i & t$Time2 == "Post" & t$Survival == "Long" & t$SampleType == "BM",]$percentage))}


[1] "Plasmablast"

Welch Two Sample t-test

data:  t[t$Var1 == i & t$Time2 == "Post" & t$Survival == "Short" & t$SampleType == "BM", ]$percentage and t[t$Var1 == i & t$Time2 == "Post" & t$Survival == "Long" & t$SampleType == "BM", ]$percentage
t = -2.3486, df = 10.989, p-value = 0.0386
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -0.0220845536 -0.0007154464
sample estimates:
  mean of x mean of y 
0.0050    0.0164 

[1] "CD8"

Welch Two Sample t-test

data:  t[t$Var1 == i & t$Time2 == "Post" & t$Survival == "Short" & t$SampleType == "BM", ]$percentage and t[t$Var1 == i & t$Time2 == "Post" & t$Survival == "Long" & t$SampleType == "BM", ]$percentage
t = -3.4217, df = 10.74, p-value = 0.005901
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -0.17482543 -0.03770791
sample estimates:
  mean of x  mean of y 
0.07133333 0.17760000 


############# Immune subtypes Welch T test comparison ######-

#Short_Pre vs Short_Post
for (i in unique(t2$Var1)) {
  print(i)
  print(t.test(t2[t2$Var1 == i & t2$Time2 == "Pre" & t2$Survival == "Short" & t2$SampleType == "BM",]$percentage, t2[t2$Var1 == i & t2$Time2 == "Post" & t2$Survival == "Short" & t2$SampleType == "BM",]$percentage))}

[1] "CD8 Naive"

Welch Two Sample t-test

data:  t2[t2$Var1 == i & t2$Time2 == "Pre" & t2$Survival == "Short" & t2$SampleType == "BM", ]$percentage and t2[t2$Var1 == i & t2$Time2 == "Post" & t2$Survival == "Short" & t2$SampleType == "BM", ]$percentage
t = 3.5, df = 2.8764, p-value = 0.04217
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  0.0003183701 0.0090149632
sample estimates:
  mean of x   mean of y 
0.007666667 0.003000000 




#Long_Pre vs Long_Post
for (i in unique(t2$Var1)) {
  print(i)
  print(t.test(t2[t2$Var1 == i & t2$Time2 == "Pre" & t2$Survival == "Long" & t2$SampleType == "BM",]$percentage, t2[t2$Var1 == i & t2$Time2 == "Post" & t2$Survival == "Long" & t2$SampleType == "BM",]$percentage))}

[1] "GMP"

Welch Two Sample t-test

data:  t2[t2$Var1 == i & t2$Time2 == "Pre" & t2$Survival == "Long" & t2$SampleType == "BM", ]$percentage and t2[t2$Var1 == i & t2$Time2 == "Post" & t2$Survival == "Long" & t2$SampleType == "BM", ]$percentage
t = 2.4407, df = 9.1118, p-value = 0.03701
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  0.00483068 0.12416932
sample estimates:
  mean of x mean of y 
0.0745    0.0100 




#Short_Pre vs Long_Pre
for (i in unique(t2$Var1)) {
  print(i)
  print(t.test(t2[t2$Var1 == i & t2$Time2 == "Pre" & t2$Survival == "Short" & t2$SampleType == "BM",]$percentage, t2[t2$Var1 == i & t2$Time2 == "Pre" & t2$Survival == "Long" & t2$SampleType == "BM",]$percentage))}



[1] "pDC"

Welch Two Sample t-test

data:  t2[t2$Var1 == i & t2$Time2 == "Pre" & t2$Survival == "Short" & t2$SampleType == "BM", ]$percentage and t2[t2$Var1 == i & t2$Time2 == "Pre" & t2$Survival == "Long" & t2$SampleType == "BM", ]$percentage
t = -3.2672, df = 10.986, p-value = 0.007514
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -0.007699272 -0.001500728
sample estimates:
  mean of x mean of y 
0.0010    0.0056 

[1] "Classical Mono"

Welch Two Sample t-test

data:  t2[t2$Var1 == i & t2$Time2 == "Pre" & t2$Survival == "Short" & t2$SampleType == "BM", ]$percentage and t2[t2$Var1 == i & t2$Time2 == "Pre" & t2$Survival == "Long" & t2$SampleType == "BM", ]$percentage
t = -2.3208, df = 9.0993, p-value = 0.04513
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -0.134499666 -0.001833667
sample estimates:
  mean of x   mean of y 
0.004333333 0.072500000 

[1] "CD4 TCM"

Welch Two Sample t-test

data:  t2[t2$Var1 == i & t2$Time2 == "Pre" & t2$Survival == "Short" & t2$SampleType == "BM", ]$percentage and t2[t2$Var1 == i & t2$Time2 == "Pre" & t2$Survival == "Long" & t2$SampleType == "BM", ]$percentage
t = 3.3212, df = 10.996, p-value = 0.00682
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  0.03970682 0.19575985
sample estimates:
  mean of x mean of y 
0.2443333 0.1266000 






#Short_Post vs Long_Post
for (i in unique(t2$Var1)) {
  print(i)
  print(t.test(t2[t2$Var1 == i & t2$Time2 == "Post" & t2$Survival == "Short" & t2$SampleType == "BM",]$percentage, t2[t2$Var1 == i & t2$Time2 == "Post" & t2$Survival == "Long" & t2$SampleType == "BM",]$percentage))}


[1] "Plasmablast"

Welch Two Sample t-test

data:  t2[t2$Var1 == i & t2$Time2 == "Post" & t2$Survival == "Short" & t2$SampleType == "BM", ]$percentage and t2[t2$Var1 == i & t2$Time2 == "Post" & t2$Survival == "Long" & t2$SampleType == "BM", ]$percentage
t = -2.3486, df = 10.989, p-value = 0.0386
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -0.0220845536 -0.0007154464
sample estimates:
  mean of x mean of y 
0.0050    0.0164 

[1] "CD8 TEM"

Welch Two Sample t-test

data:  t2[t2$Var1 == i & t2$Time2 == "Post" & t2$Survival == "Short" & t2$SampleType == "BM", ]$percentage and t2[t2$Var1 == i & t2$Time2 == "Post" & t2$Survival == "Long" & t2$SampleType == "BM", ]$percentage
t = -3.1698, df = 10.355, p-value = 0.009581
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -0.07444534 -0.01315466
sample estimates:
  mean of x mean of y 
0.0180    0.0618 

[1] "CD8 TEMRA"

Welch Two Sample t-test

data:  t2[t2$Var1 == i & t2$Time2 == "Post" & t2$Survival == "Short" & t2$SampleType == "BM", ]$percentage and t2[t2$Var1 == i & t2$Time2 == "Post" & t2$Survival == "Long" & t2$SampleType == "BM", ]$percentage
t = -2.9185, df = 10.425, p-value = 0.01472
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -0.09681783 -0.01324883
sample estimates:
  mean of x  mean of y 
0.03766667 0.09270000 



[1] "CD56bright NK"

Welch Two Sample t-test

data:  t2[t2$Var1 == i & t2$Time2 == "Post" & t2$Survival == "Short" & t2$SampleType == "BM", ]$percentage and t2[t2$Var1 == i & t2$Time2 == "Post" & t2$Survival == "Long" & t2$SampleType == "BM", ]$percentage
t = -2.3488, df = 10.364, p-value = 0.03989
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -0.0152939557 -0.0004393776
sample estimates:
  mean of x   mean of y 
0.003333333 0.011200000 



### PLOT Boxplot GRAPH OF IMMUNE FREQUENCY####
t$Patient <- factor(t$Patient, levels = c("24780","122759","113446","122739","121535","121474","119837","121015","122976","114152","119325","118877","115585"))
PatientColor <- c("#111111","#555555","#999999","#f1ddf1","#eaccea","#e3bbe3","#ddaadd","#c699c6","#de94f5","#cb53ef","#be29ec","#9820bc","#851ca5")
PatientColorBox <- c("#999999","#cb53ef","#851ca5","#9820bc","#be29ec","#c699c6","#ddaadd","#e3bbe3","#eaccea","#f1ddf1","#555555","#de94f5","#111111")
t$Time2 <- factor(t$Time2, levels = c("Pre","C1D8","Post","C2D8"))
t$Survival <- factor(t$Survival, levels = c("Short","Long"))
t$TimeSurvival <- paste(t$Survival,t$Time2,sep = "_")
t$TimeSurvival <- factor(t$TimeSurvival, levels = c("Short_Pre","Short_C1D8","Short_Post","Short_C2D8","Long_Pre","Long_C1D8","Long_Post"))


p<-ggplot(t[t$SampleType == "BM",])+geom_boxplot(aes(Time2,percentage, group = interaction(Time2, Survival), color=Survival),outlier.shape = NA,width = 0.4)+geom_line(aes(Time2,percentage, group = Patient, color=Patient),lwd=1,alpha=0.3)+geom_point(aes(Time2,percentage,group=interaction(Time2, Survival), color=factor(Patient)),size=3,position = position_dodge(width=0.4))+ggtitle("Immune Frequency Dynamics (Broad)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Treatment", y = "Frequency")+scale_fill_manual(values =c("#555555","#c699c6",PatientColor))+scale_color_manual(values =c("#555555","#c699c6",PatientColor))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
p+ facet_wrap( ~ Var1, scales="free",ncol=8)

ggplot(t[t$SampleType == "BM",])+geom_boxplot(aes(Var1,percentage, group = interaction(Var1, TimeSurvival), color=TimeSurvival,fill=TimeSurvival,alpha=0.3),outlier.shape = NA,width = 0.8)+geom_point(aes(Var1,percentage, group = interaction(Var1, TimeSurvival), color=TimeSurvival),size=3,position = position_dodge(width=0.8))+ggtitle("Immune Frequency Dynamics (Broad)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Immune Subtype", y = "Frequency")+scale_fill_manual(values =c("#999999","#111111","#ddaadd","#851ca5"))+scale_color_manual(values =c("#999999","#111111","#ddaadd","#851ca5"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))






##Anno_Specific
t2$Patient <- factor(t2$Patient, levels = c("24780","122759","113446","122739","121535","121474","119837","121015","122976","114152","119325","118877","115585"))
PatientColor <- c("#111111","#555555","#999999","#f1ddf1","#eaccea","#e3bbe3","#ddaadd","#c699c6","#de94f5","#cb53ef","#be29ec","#9820bc","#851ca5")
PatientColorBox <- c("#999999","#cb53ef","#851ca5","#9820bc","#be29ec","#c699c6","#ddaadd","#e3bbe3","#eaccea","#f1ddf1","#555555","#de94f5","#111111")
t2$Time2 <- factor(t2$Time2, levels = c("Pre","C1D8","Post","C2D8"))
t2$Survival <- factor(t2$Survival, levels = c("Short","Long"))
t2$TimeSurvival <- paste(t2$Survival,t2$Time2,sep = "_")
t2$TimeSurvival <- factor(t2$TimeSurvival, levels = c("Short_Pre","Short_C1D8","Short_Post","Short_C2D8","Long_Pre","Long_C1D8","Long_Post"))
t2 <- t2[order(t2$Patient),]

p<-ggplot(t2[t2$SampleType == "BM",])+geom_boxplot(aes(Time2,percentage, group = interaction(Time2, Survival), color=Survival),outlier.shape = NA,width = 0.4)+geom_line(aes(Time2,percentage, group = Patient, color=Patient),lwd=1,alpha=0.3)+geom_point(aes(Time2,percentage,group=interaction(Time2, Survival), color=factor(Patient)),size=3,position = position_dodge(width=0.4))+ggtitle("Immune Frequency Dynamics (Specific)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Treatment", y = "Frequency")+scale_fill_manual(values =c("#555555","#c699c6",PatientColor))+scale_color_manual(values =c("#555555","#c699c6",PatientColor))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
p+ facet_wrap( ~ Var1, scales="free",ncol=8)

ggplot(t2[t2$SampleType == "BM",])+geom_boxplot(aes(Var1,percentage, group = interaction(Var1, TimeSurvival), color=TimeSurvival,fill=TimeSurvival,alpha=0.3),outlier.shape = NA,width = 0.8)+geom_point(aes(Var1,percentage, group = interaction(Var1, TimeSurvival), color=TimeSurvival),size=2,position = position_dodge(width=0.8))+ggtitle("Immune Frequency Dynamics (Specific)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Immune Subtype", y = "Frequency")+scale_fill_manual(values =c("#999999","#111111","#ddaadd","#851ca5"))+scale_color_manual(values =c("#999999","#111111","#ddaadd","#851ca5"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))

ggplot(t2[t2$SampleType == "BM",])+geom_boxplot(aes(Var1,percentage, group = interaction(Var1, TimeSurvival), color=TimeSurvival,fill=TimeSurvival,alpha=0.3),outlier.shape = NA,width = 0.8)+geom_point(aes(Var1,percentage, group = interaction(Var1, TimeSurvival), color=Patient),size=2,position = position_dodge(width=0.8))+ggtitle("Immune Count Dynamics (Specifc)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Immune Subtype", y = "Number of Cells")+scale_fill_manual(values =c("#999999","#111111","#ddaadd","#851ca5",mypalette,mypalette2))+scale_color_manual(values =c("#999999","#111111","#ddaadd","#851ca5",mypalette,mypalette2))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))

ggplot(t2[t2$SampleType == "BM",])+geom_boxplot(aes(Var1,Freq, group = interaction(Var1, TimeSurvival), color=TimeSurvival,fill=TimeSurvival,alpha=0.3),outlier.shape = NA,width = 0.8)+geom_point(aes(Var1,Freq, group = interaction(Var1, TimeSurvival), color=Patient),size=2,position = position_dodge(width=0.8))+ggtitle("Immune Count Dynamics (Specifc)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Immune Subtype", y = "Number of Cells")+scale_fill_manual(values =c("#999999","#111111","#ddaadd","#851ca5",mypalette,mypalette2))+scale_color_manual(values =c("#999999","#111111","#ddaadd","#851ca5",mypalette,mypalette2))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+geom_hline(yintercept = 20, lwd = 2)

ggplot(t2[t2$SampleType == "BM" & t2$Var1 %in% c("Intermediate Mono","Non-classical Mono","mDC","pDC","Neutrophil","Megakaryocyte","Stromal"),])+geom_boxplot(aes(Var1,percentage, group = interaction(Var1, TimeSurvival), color=TimeSurvival,fill=TimeSurvival,alpha=0.3),outlier.shape = NA,width = 0.8)+geom_point(aes(Var1,percentage, group = interaction(Var1, TimeSurvival), color=TimeSurvival),size=2,position = position_dodge(width=0.8))+ggtitle("Immune Frequency Dynamics (Broad)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Immune Subtype", y = "Frequency")+coord_cartesian(ylim=c(0,0.08))+scale_fill_manual(values =c("#999999","#111111","#ddaadd","#851ca5"))+scale_color_manual(values =c("#999999","#111111","#ddaadd","#851ca5"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))



##################### USE WILCOX RANK SUM TEST FOR DEG ANALYSIS EXPLORATION #######################
BM.combined.sct <- LoadSeuratRds("BM.combined.sct_2025.rds") ######## LOAD ONLY THE BM SEURAT DATASET TO REDUCE MEMORY CONSUMPTION AND SERVER CRASH

BM.combined.sct$SampleInfo_Type_Time <- paste(BM.combined.sct$Patient,BM.combined.sct$SampleType,BM.combined.sct$Time,sep="-")
BM.combined.sct$SampleInfo_Type_Time<- factor(BM.combined.sct$SampleInfo_Type_Time, levels = c("24780-BM-Scr","24780-BM-PC2","122759-BM-Scr","122759-BM-PC2","113446-BM-Scr","113446-BM-PC2","122739-BM-Scr","122739-BM-PC2","121535-BM-Scr","121535-BM-PC2","121474-BM-Scr","121474-BM-PC2","119837-BM-Scr","119837-BM-PC2","121015-BM-Scr","121015-BM-PC2","122976-BM-Scr","122976-BM-PC2","114152-BM-Scr","114152-BM-PC2","119325-BM-Scr","119325-BM-PC2","118877-BM-Scr","118877-BM-PC2","115585-BM-Scr","115585-BM-PC2","24780-PBMC-Scr","24780-PBMC-C2D1","24780-PBMC-C2D8","121474-PBMC-Scr","121474-PBMC-C1D8","121474-PBMC-C2D1","121015-PBMC-Scr","121015-PBMC-C1D8","121015-PBMC-C2D1","118877-PBMC-Scr","118877-PBMC-C1D8","118877-PBMC-C2D1"))

BM.combined.sct@meta.data$PatientTime <- paste(BM.combined.sct@meta.data$Patient,BM.combined.sct@meta.data$Time2, sep="_")
BM.combined.sct@meta.data$TimeSurvival <- paste(BM.combined.sct@meta.data$Survival,BM.combined.sct@meta.data$Time2,sep = "_")

BM.combined.sct.BM <- subset(BM.combined.sct, subset = SampleType == "BM")
BM.combined.sct.PBMC <- subset(BM.combined.sct, subset = SampleType == "PBMC")


DefaultAssay(BM.combined.sct.BM) = "SCT"
Idents(BM.combined.sct.BM) <- "manualAnno_specific"
BM.combined.sct.BM$BroadAnno_TimeSurvival <- paste(BM.combined.sct.BM$manualAnno_specific,BM.combined.sct.BM$TimeSurvival,sep = "-")
BM.combined.sct.BM$SpecificAnno_TimeSurvival <- paste(BM.combined.sct.BM$manualAnno_specific2,BM.combined.sct.BM$TimeSurvival,sep = "-")


DEGs_SpecificAnno_wilcox = data.frame()
for (i in as.character(unique(BM.combined.sct.BM$manualAnno_specific2))) 
{tryCatch({
  Idents(BM.combined.sct.BM) <- "SpecificAnno_TimeSurvival"
  print(i)
  DEGS <- FindMarkers(BM.combined.sct.BM, assay = "SCT", ident.1 = paste(i,"Short_Pre",sep="-"), ident.2 = paste(i,"Short_Post",sep="-"),
                      verbose = FALSE,recorrect_umi=FALSE,logfc.threshold = 0, min.pct = 0, test.use="wilcox")  
  DEGS$manualAnno_specific2 <- i
  DEGS$comparison <- "Short-Pre_vs_Post"
  DEGS$Gene <- rownames(DEGS)
  DEGs_SpecificAnno_wilcox <- rbind(DEGs_SpecificAnno_wilcox,DEGS)
  DEGS <- FindMarkers(BM.combined.sct.BM, assay = "SCT", ident.1 = paste(i,"Long_Pre",sep="-"), ident.2 = paste(i,"Long_Post",sep="-"),verbose = FALSE,recorrect_umi=FALSE,logfc.threshold = 0, min.pct = 0, test.use="wilcox")   
  DEGS$manualAnno_specific2 <- i
  DEGS$comparison <- "Long-Pre_vs_Post"
  DEGS$Gene <- rownames(DEGS)
  DEGs_SpecificAnno_wilcox <- rbind(DEGs_SpecificAnno_wilcox,DEGS)
  DEGS <- FindMarkers(BM.combined.sct.BM, assay = "SCT", ident.1 = paste(i,"Short_Pre",sep="-"), ident.2 = paste(i,"Long_Pre",sep="-"),verbose = FALSE,recorrect_umi=FALSE,logfc.threshold = 0, min.pct = 0, test.use="wilcox")  
  DEGS$manualAnno_specific2 <- i
  DEGS$comparison <- "Short_vs_Long-Pre"
  DEGS$Gene <- rownames(DEGS)
  DEGs_SpecificAnno_wilcox <- rbind(DEGs_SpecificAnno_wilcox,DEGS)
  DEGS <- FindMarkers(BM.combined.sct.BM, assay = "SCT", ident.1 = paste(i,"Short_Post",sep="-"), ident.2 = paste(i,"Long_Post",sep="-"),verbose = FALSE,recorrect_umi=FALSE,logfc.threshold = 0, min.pct = 0, test.use="wilcox")  
  DEGS$manualAnno_specific2 <- i
  DEGS$comparison <- "Short_vs_Long-Post"
  DEGS$Gene <- rownames(DEGS)
  DEGs_SpecificAnno_wilcox <- rbind(DEGs_SpecificAnno_wilcox,DEGS)
}, 
error=function(e){cat("ERROR :",conditionMessage(e), "\n")})}


DEGs_SpecificAnno_wilcox[DEGs_SpecificAnno_wilcox$p_val_adj == 0,]$p_val_adj <- .Machine$double.xmin ### Change adjust pvalue to lowest number possible in R (to avoid Inf when doing log(adjP))

DEGs_SpecificAnno.WL <- DEGs_SpecificAnno_wilcox[-grep("MT-",rownames(DEGs_SpecificAnno_wilcox)),]#remove MT genes from DEGs
DEGs_SpecificAnno.WL <- na.omit(DEGs_SpecificAnno.WL)
DEGs_SpecificAnno.WL <- DEGs_SpecificAnno.WL[DEGs_SpecificAnno.WL$p_val_adj < 0.01,]
DEGs_SpecificAnno.WL <- DEGs_SpecificAnno.WL[DEGs_SpecificAnno.WL$pct.1 > 0.20 | DEGs_SpecificAnno.WL$pct.2 > 0.20,] #filter genes for expression captured in at least 20% of cells.


for (i in unique(DEGs_SpecificAnno.WL$comparison)) {
  print (i)
  print (table(DEGs_SpecificAnno.WL[DEGs_SpecificAnno.WL$comparison ==i,]$manualAnno_specific2))
}



##### PERFORM DIFFERENTIAL ADT ABUNDANCE
DefaultAssay(BM.combined.sct.BM) = "ADT2"

dADTs_SpecificAnno_wilcox = data.frame()
for (i in as.character(unique(BM.combined.sct.BM$manualAnno_specific2))) 
{tryCatch({
  Idents(BM.combined.sct.BM) <- "SpecificAnno_TimeSurvival"
  print(i)
  dADTS <- FindMarkers(BM.combined.sct.BM, assay = "ADT2", ident.1 = paste(i,"Short_Pre",sep="-"), ident.2 = paste(i,"Short_Post",sep="-"),
                       verbose = FALSE,recorrect_umi=FALSE,logfc.threshold = 0, min.pct = 0, test.use="wilcox")  
  dADTS$manualAnno_specific2 <- i
  dADTS$comparison <- "Short-Pre_vs_Post"
  dADTS$Gene <- rownames(dADTS)
  dADTs_SpecificAnno_wilcox <- rbind(dADTs_SpecificAnno_wilcox,dADTS)
  dADTS <- FindMarkers(BM.combined.sct.BM, assay = "ADT2", ident.1 = paste(i,"Long_Pre",sep="-"), ident.2 = paste(i,"Long_Post",sep="-"),verbose = FALSE,recorrect_umi=FALSE,logfc.threshold = 0, min.pct = 0, test.use="wilcox")   
  dADTS$manualAnno_specific2 <- i
  dADTS$comparison <- "Long-Pre_vs_Post"
  dADTS$Gene <- rownames(dADTS)
  dADTs_SpecificAnno_wilcox <- rbind(dADTs_SpecificAnno_wilcox,dADTS)
  dADTS <- FindMarkers(BM.combined.sct.BM, assay = "ADT2", ident.1 = paste(i,"Short_Pre",sep="-"), ident.2 = paste(i,"Long_Pre",sep="-"),verbose = FALSE,recorrect_umi=FALSE,logfc.threshold = 0, min.pct = 0, test.use="wilcox")  
  dADTS$manualAnno_specific2 <- i
  dADTS$comparison <- "Short_vs_Long-Pre"
  dADTS$Gene <- rownames(dADTS)
  dADTs_SpecificAnno_wilcox <- rbind(dADTs_SpecificAnno_wilcox,dADTS)
  dADTS <- FindMarkers(BM.combined.sct.BM, assay = "ADT2", ident.1 = paste(i,"Short_Post",sep="-"), ident.2 = paste(i,"Long_Post",sep="-"),verbose = FALSE,recorrect_umi=FALSE,logfc.threshold = 0, min.pct = 0, test.use="wilcox")  
  dADTS$manualAnno_specific2 <- i
  dADTS$comparison <- "Short_vs_Long-Post"
  dADTS$Gene <- rownames(dADTS)
  dADTs_SpecificAnno_wilcox <- rbind(dADTs_SpecificAnno_wilcox,dADTS)
}, 
error=function(e){cat("ERROR :",conditionMessage(e), "\n")})}



dADTs_SpecificAnno.WL <- na.omit(dADTs_SpecificAnno_wilcox)
dADTs_SpecificAnno.WL <- dADTs_SpecificAnno.WL[dADTs_SpecificAnno.WL$p_val_adj < 0.01,]
dADTs_SpecificAnno.WL <- dADTs_SpecificAnno.WL[dADTs_SpecificAnno.WL$pct.1 > 0.20 | dADTs_SpecificAnno.WL$pct.2 > 0.20,] #filter genes for expression captured in at least 20% of cells.

for (i in unique(dADTs_SpecificAnno.WL$comparison)) {
  print (i)
  print (table(dADTs_SpecificAnno.WL[dADTs_SpecificAnno.WL$comparison ==i,]$manualAnno_specific2))
}



### MAKE HEATMAP OF SIGNIFICANT LOG FOLD CHANGES OF IMMUNE CHEKCPOINT GENES ####
ICI <- c("CD279","CD274-B7H1-PDL1","CD223-LAG3","CD272-BTLA","CD366-Tim3","TIGIT-VSTM3","CD96-TACTILE","CD155-PVR")

dADTs_SpecificAnno.WL_ICI <- dADTs_SpecificAnno.WL[dADTs_SpecificAnno.WL$Gene %in% ICI,]
dADTs_SpecificAnno.WL_ICI <- dADTs_SpecificAnno.WL_ICI[dADTs_SpecificAnno.WL_ICI$manualAnno_specific2 %in% c("EPC","MEP","HSC","LMPP","GMP","Classical Mono","Intermediate Mono","Non-classical Mono","mDC","pDC","Naive B","Memory B","Plasmablast","Vd1 gdT","Vd2 gdT","CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL","CD4 Treg","CD8 Naive","CD8 TCM","CD8 TEM","CD8 TEMRA","MAIT","CD56dim NK","CD56bright NK"),]
                                                      
dADTs_SpecificAnno.WL_ICI$Gene <- factor(dADTs_SpecificAnno.WL_ICI$Gene, levels = c(ICI))
dADTs_SpecificAnno.WL_ICI$manualAnno_specific2 <- factor(dADTs_SpecificAnno.WL_ICI$manualAnno_specific2, levels = c("EPC","MEP","HSC","LMPP","GMP","Classical Mono","Intermediate Mono","Non-classical Mono","mDC","pDC","Naive B","Memory B","Plasmablast","Vd1 gdT","Vd2 gdT","CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL","CD4 Treg","CD8 Naive","CD8 TCM","CD8 TEM","CD8 TEMRA","MAIT","CD56dim NK","CD56bright NK"))
                                                       
dADTs_SpecificAnno.WL_ICI[dADTs_SpecificAnno.WL_ICI$avg_log2FC > 2,]$avg_log2FC <- 2 #for visualization purposes, change anything > 2 to 2
dADTs_SpecificAnno.WL_ICI[dADTs_SpecificAnno.WL_ICI$avg_log2FC < -2,]$avg_log2FC <- -2 #for visualization purposes, change anything < -2 to -2
                                                       
ggplot(dADTs_SpecificAnno.WL_ICI[dADTs_SpecificAnno.WL_ICI$comparison == "Short_vs_Long-Pre",], aes(y = manualAnno_specific2, x = Gene, fill = avg_log2FC))+geom_tile()+scale_y_discrete(limits=rev)+ scale_fill_gradient2(low = "#8247a3", high = "#5f5e60",limits = c(-2,2)) +theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ geom_point(aes(size =-log10(p_val_adj)))
                                                       

ggplot(dADTs_SpecificAnno.WL_ICI[dADTs_SpecificAnno.WL_ICI$comparison == "Short_vs_Long-Pre"& dADTs_SpecificAnno.WL_ICI$Gene %in% c("CD279","CD274-B7H1-PDL1","CD223-LAG3","CD366-Tim3","CD272-BTLA","TIGIT-VSTM3","CD96-TACTILE","CD155-PVR"),], aes(y = manualAnno_specific2, x = Gene, fill = avg_log2FC))+geom_tile()+scale_y_discrete(limits=rev)+ scale_fill_gradient2(low = "#8247a3", high = "#5f5e60",limits = c(-2,2)) +theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ geom_point(aes(size =-log10(p_val_adj)))





######### MAKE VOLCANO PLOTS OF WILCOX DEGS ###############
volcdata <- DEGs_SpecificAnno_wilcox[DEGs_SpecificAnno_wilcox$manualAnno_specific2 == "MEP" & DEGs_SpecificAnno_wilcox$comparison == "Short_vs_Long-Pre" & (DEGs_SpecificAnno_wilcox$pct.1 > 0.15 | DEGs_SpecificAnno_wilcox$pct.2 > 0.15),]
EnhancedVolcano(volcdata,
                lab = volcdata$Gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                #selectLab = c("ANXA1","C12orf75",'CFD','CP','AC044893.1',"CLK3","PROK2","CD96", "CDKN2C", "CDKN2D","MPO","CPT1C","RGPD2","TYROBP","LINC00355","TCEAL2"), #Pre
                #selectLab = c("TYROBP","MGST1","CFD","CD96","CLK3","MPO","EGFL6","FBLN5","PRSS2","CDKN2C","CDKN2D","CP","LGALS1","CA2","LMNA","AZU1","PROK2","MFAP4","AC244502.1","RAB7B"), #post
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'P.adj'),
                pCutoff = 0.01,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 5.0,
                colAlpha = 0.8,
                labFace = 'bold',
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                col=c('grey', 'grey', 'grey', 'red3'),
                colConnectors = 'black',
                drawConnectors = F,
                widthConnectors = 0.75)



volcdata <- DEGs_SpecificAnno_wilcox[DEGs_SpecificAnno_wilcox$manualAnno_specific2 == "HSC" & DEGs_SpecificAnno_wilcox$comparison == "Short_vs_Long-Pre" & (DEGs_SpecificAnno_wilcox$pct.1 > 0.15 | DEGs_SpecificAnno_wilcox$pct.2 > 0.15),]
EnhancedVolcano(volcdata,
                lab = volcdata$Gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                selectLab = c("ANXA1","C12orf75",'CFD','CP','AC044893.1',"CLK3","PROK2","CD96", "CDKN2C", "CDKN2D","MPO","CPT1C","RGPD2","TYROBP","LINC00355","TCEAL2"), #Pre
                #selectLab = c("TYROBP","MGST1","CFD","CD96","CLK3","MPO","EGFL6","FBLN5","PRSS2","CDKN2C","CDKN2D","CP","LGALS1","CA2","LMNA","AZU1","PROK2","MFAP4","AC244502.1","RAB7B"), #post
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'P.adj'),
                pCutoff = 0.01,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 5.0,
                colAlpha = 0.8,
                labFace = 'bold',
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                col=c('grey', 'grey', 'grey', 'red3'),
                colConnectors = 'black',
                drawConnectors = F,
                widthConnectors = 0.75)


volcdata <- DEGs_SpecificAnno_wilcox[DEGs_SpecificAnno_wilcox$manualAnno_specific2 == "LMPP" & DEGs_SpecificAnno_wilcox$comparison == "Short_vs_Long-Pre"& (DEGs_SpecificAnno_wilcox$pct.1 > 0.15 | DEGs_SpecificAnno_wilcox$pct.2 > 0.15),]

EnhancedVolcano(volcdata,
                lab = volcdata$Gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                selectLab = c("PTPRCAP","TXNIP","ADA","YWHAB","GNAI1","CPXM1","GSDMD","GNAI1","CA2", "TRBC1","CLK3","CD96","CDKN2C","CXCL2","CXCL3","GPC6","CXCL8","PCAT18","STAB1"),#pre
                #selectLab = c("SELL","SMIM24","GYPC","IGFBP7","FABP5","SHMT2","CD96","CXCL2","NKX2-5","S100A16","CST3","KLF2","CD68","RHOB","BGLAP"),#Post
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'P.adj'),
                pCutoff = 0.01,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 5.0,
                colAlpha = 0.8,
                labFace = 'bold',
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                col=c('grey', 'grey', 'grey', 'red3'),
                colConnectors = 'black',
                drawConnectors = F,
                widthConnectors = 0.75)

volcdata <- DEGs_SpecificAnno_wilcox[DEGs_SpecificAnno_wilcox$manualAnno_specific2 == "GMP" & DEGs_SpecificAnno_wilcox$comparison == "Short_vs_Long-Pre"& (DEGs_SpecificAnno_wilcox$pct.1 > 0.15 | DEGs_SpecificAnno_wilcox$pct.2 > 0.15),]

EnhancedVolcano(volcdata,
                lab = volcdata$Gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                selectLab = c("MNDA","TXNIP","S100A8","S100A9","GTSF1","S100A12","CA2","CLK3","CD96","CDKN2C","CDKN2D","CXCL2","CXCL3","CCL5","PRSS21","LGALS3BP","CD69","AL392086.3", "IL7","SPINK2","KCNE5","BEX2"),#pre
                #selectLab = c("S100A8","S100A9","COTL1","TSPO","LYZ","CD96","CA2","CLK3","CDKN2C","CXCL3","SH2D1A","LGALS3BP","MSI2","SPINK2","MNDA"),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'P.adj'),
                pCutoff = 0.01,
                FCcutoff = 1.0,
                pointSize = 2.0,
                labSize = 5.0,
                colAlpha = 0.8,
                labFace = 'bold',
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                col=c('grey', 'grey', 'grey', 'red3'),
                colConnectors = 'black',
                drawConnectors = F,
                widthConnectors = 0.75)


volcdata <- DEGs_SpecificAnno_wilcox[DEGs_SpecificAnno_wilcox$manualAnno_specific2 == "mDC" & DEGs_SpecificAnno_wilcox$comparison == "Short_vs_Long-Pre"& (DEGs_SpecificAnno_wilcox$pct.1 > 0.15 | DEGs_SpecificAnno_wilcox$pct.2 > 0.15),]

EnhancedVolcano(volcdata,
                lab = volcdata$Gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                #selectLab = c("RPS4Y1", "HLA-DRB5","CEBPD","CD48","HLA-DQA2","CD1C","CD1D","CD1E","LMNA", "NKX2-5", "S100A16","BASP1"),
                selectLab = c("RPS4Y1", "HLA-DRB5","CEBPD","CD48","HLA-DQA2","CD1C","CD1D","CD1E","LMNA", "NKX2-5", "CFD","RHOB"),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'P.adj'),
                pCutoff = 0.0001,
                FCcutoff = 1.0,
                pointSize = 2.0,
                labSize = 5.0,
                colAlpha = 0.8,
                labFace = 'bold',
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                col=c('grey', 'grey', 'grey', 'red3'),
                colConnectors = 'black',
                drawConnectors = F,
                widthConnectors = 0.75)


volcdata <- DEGs_SpecificAnno_wilcox[DEGs_SpecificAnno_wilcox$manualAnno_specific2 == "Classical Mono" & DEGs_SpecificAnno_wilcox$comparison == "Short_vs_Long-Pre"& (DEGs_SpecificAnno_wilcox$pct.1 > 0.15 | DEGs_SpecificAnno_wilcox$pct.2 > 0.15),]
volcdata <- DEGs_SpecificAnno_wilcox[DEGs_SpecificAnno_wilcox$manualAnno_specific2 == "Classical Mono" & DEGs_SpecificAnno_wilcox$comparison == "Short_vs_Long-Post"& (DEGs_SpecificAnno_wilcox$pct.1 > 0.15 | DEGs_SpecificAnno_wilcox$pct.2 > 0.15),]

EnhancedVolcano(volcdata,
                lab = volcdata$Gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                #selectLab = c("RPS4Y1", "HLA-DRB5","CEBPD","CD48","HLA-DQA2","CD1C","CD1D","CD1E","LMNA", "NKX2-5", "S100A16","BASP1"),
                #selectLab = c("RPS4Y1", "HLA-DRB5","CEBPD","CD48","HLA-DQA2","CD1C","CD1D","CD1E","LMNA", "NKX2-5", "CFD","RHOB"),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'P.adj'),
                pCutoff = 0.0001,
                FCcutoff = 1.0,
                pointSize = 2.0,
                labSize = 5.0,
                colAlpha = 0.8,
                labFace = 'bold',
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                col=c('grey', 'grey', 'grey', 'red3'),
                colConnectors = 'black',
                drawConnectors = F,
                widthConnectors = 0.75)





### Look for reoccurring DEGs in various subsets 
#group 1 - HSC, LMPP, GMP
#group 2 - GMP, Classical, Intermediate, Non-classical mono, mDC, pDC
#group 3 - CD4 subsets
#group 4 - CD8 subsets
#group 5 - NK, ydT
#group 6 - B cells

### Group 1 
DEGs_SpecificAnno.WL1 <- DEGs_SpecificAnno.WL[DEGs_SpecificAnno.WL$manualAnno_specific2 %in% c("GMP","LMPP","HSC") & DEGs_SpecificAnno.WL$comparison %in% c("Short_vs_Long-Pre"),]
dt <- data.frame(table(DEGs_SpecificAnno.WL1$Gene))
DEGs_SpecificAnno.WL2 <- DEGs_SpecificAnno.WL1[order(DEGs_SpecificAnno.WL1$Gene),]
DEGs_SpecificAnno.WL2<- DEGs_SpecificAnno.WL2[DEGs_SpecificAnno.WL2$Gene %in% dt[dt$Freq >= 2,]$Var1 & abs(DEGs_SpecificAnno.WL2$avg_log2FC) > 1,]
DEGs_SpecificAnno.WL1[DEGs_SpecificAnno.WL1$Gene %in% DEGs_SpecificAnno.WL2$Gene,]

### CD8
DEGs_SpecificAnno.WL1 <- DEGs_SpecificAnno.WL[DEGs_SpecificAnno.WL$manualAnno_specific2 %in% c("CD8 Naive","CD8 TCM","CD8 TEM", "CD8 TEMRA") & DEGs_SpecificAnno.WL$comparison %in% c("Short_vs_Long-Post"),]
dt <- data.frame(table(DEGs_SpecificAnno.WL1$Gene))
DEGs_SpecificAnno.WL2 <- DEGs_SpecificAnno.WL1[order(DEGs_SpecificAnno.WL1$Gene),]
DEGs_SpecificAnno.WL2<- DEGs_SpecificAnno.WL1[DEGs_SpecificAnno.WL1$Gene %in% dt[dt$Freq >= 1,]$Var1 & abs(DEGs_SpecificAnno.WL1$avg_log2FC) > 1,]
DEGs_SpecificAnno.WL2[order(DEGs_SpecificAnno.WL2$avg_log2FC),]

#CD4 Naive, TCM
DEGs_SpecificAnno.WL1 <- DEGs_SpecificAnno.WL[DEGs_SpecificAnno.WL$manualAnno_specific2 %in% c("CD4 Naive","CD4 TCM") & DEGs_SpecificAnno.WL$comparison %in% c("Short_vs_Long-Pre"),]
dt <- data.frame(table(DEGs_SpecificAnno.WL1$Gene))
DEGs_SpecificAnno.WL2 <- DEGs_SpecificAnno.WL1[order(DEGs_SpecificAnno.WL1$Gene),]
DEGs_SpecificAnno.WL2<- DEGs_SpecificAnno.WL1[DEGs_SpecificAnno.WL1$Gene %in% dt[dt$Freq >= 1,]$Var1 & abs(DEGs_SpecificAnno.WL1$avg_log2FC) > 1,]
DEGs_SpecificAnno.WL2[order(DEGs_SpecificAnno.WL2$avg_log2FC),]






############# PERFORM PSEUDOBULK by PATIENTS AND PERFORM DEGs with DESeq2 (https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html, https://www.nature.com/articles/s41467-021-25960-2) #########################

#START NEW SESSION TO RESET MEMORY CONSUMPTION FROM previous DEG calling
BM.combined.sct <- LoadSeuratRds("BM.combined.sct_2025.rds") ######## LOAD ONLY THE BM SEURAT DATASET TO REDUCE MEMORY CONSUMPTION AND SERVER CRASH

BM.combined.sct$SampleInfo_Type_Time <- paste(BM.combined.sct$Patient,BM.combined.sct$SampleType,BM.combined.sct$Time,sep="-")
BM.combined.sct$SampleInfo_Type_Time<- factor(BM.combined.sct$SampleInfo_Type_Time, levels = c("24780-BM-Scr","24780-BM-PC2","122759-BM-Scr","122759-BM-PC2","113446-BM-Scr","113446-BM-PC2","122739-BM-Scr","122739-BM-PC2","121535-BM-Scr","121535-BM-PC2","121474-BM-Scr","121474-BM-PC2","119837-BM-Scr","119837-BM-PC2","121015-BM-Scr","121015-BM-PC2","122976-BM-Scr","122976-BM-PC2","114152-BM-Scr","114152-BM-PC2","119325-BM-Scr","119325-BM-PC2","118877-BM-Scr","118877-BM-PC2","115585-BM-Scr","115585-BM-PC2","24780-PBMC-Scr","24780-PBMC-C2D1","24780-PBMC-C2D8","121474-PBMC-Scr","121474-PBMC-C1D8","121474-PBMC-C2D1","121015-PBMC-Scr","121015-PBMC-C1D8","121015-PBMC-C2D1","118877-PBMC-Scr","118877-PBMC-C1D8","118877-PBMC-C2D1"))

BM.combined.sct@meta.data$PatientTime <- paste(BM.combined.sct@meta.data$Patient,BM.combined.sct@meta.data$Time2, sep="_")
BM.combined.sct@meta.data$TimeSurvival <- paste(BM.combined.sct@meta.data$Survival,BM.combined.sct@meta.data$Time2,sep = "_")

BM.combined.sct.BM <- subset(BM.combined.sct, subset = SampleType == "BM")
BM.combined.sct.PBMC <- subset(BM.combined.sct, subset = SampleType == "PBMC")

##### BONE MARROW PSEUDOBULKING ######
DefaultAssay(BM.combined.sct.BM) = "RNA"
BM.combined.sct.BM <- JoinLayers(BM.combined.sct.BM)

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- BM.combined.sct.BM[["RNA"]]$counts 

# Set up metadata as desired for aggregation and DE analysis
BM.combined.sct.BM$sample_id <- paste(BM.combined.sct.BM$Patient,BM.combined.sct.BM$TimeSurvival, sep = "_")

BM.combined.sct.BM$sample_id <- factor(BM.combined.sct.BM$sample_id)


metadata <- BM.combined.sct.BM@meta.data


# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)


# Single-cell RNA-seq analysis - Pseudobulk DE analysis with DESeq2library(tidyverse)
library(cowplot)
library(edgeR)
library(Matrix)
library(Matrix.utils) #remotes::install_github("cvarrichio/Matrix.utils")
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)
library(tibble)

# Extract unique names of clusters (= levels of cluster_id factor variable)
cluster_names_broad <- unique(colData(sce)$manualAnno_specific)
cluster_names_broad 

cluster_names_specific <- unique(colData(sce)$manualAnno_specific2)
cluster_names_specific

# Extract unique names of samples (= levels of sample_id factor variable)
sample_names <- levels(colData(sce)$sample_id)
sample_names


# Subset metadata to include only the variables you want to aggregate across (here, we want to aggregate by sample and by cluster)
groups_broad <- colData(sce)[, c("manualAnno_specific","sample_id")]
head(groups_broad)

groups_specific <- colData(sce)[, c("manualAnno_specific2","sample_id")]
head(groups_specific)


# Aggregate across cluster-sample groups
# transposing row/columns to have cell_ids as row names matching those of groups
aggr_counts_broad <- aggregate.Matrix(t(counts(sce)), 
                                      groupings = groups_broad, fun = "sum") 
aggr_counts_specific <- aggregate.Matrix(t(counts(sce)), 
                                         groupings = groups_specific, fun = "sum") 


# Transpose aggregated matrix to have genes as rows and samples as columns
aggr_counts_broad <- t(aggr_counts_broad)
aggr_counts_specific <- t(aggr_counts_specific)



# Loop over all cell types to extract corresponding counts, and store information in a list

## Initiate empty list
counts_ls_broad <- list()

for (i in 1:length(cluster_names_broad)) {
  
  ## Extract indexes of columns in the global matrix that match a given cluster
  column_idx <- which(tstrsplit(colnames(aggr_counts_broad), "_")[[1]] == cluster_names_broad[i])
  
  ## Store corresponding sub-matrix as one element of a list
  counts_ls_broad[[i]] <- aggr_counts_broad[, column_idx]
  names(counts_ls_broad)[i] <- as.character(cluster_names_broad[i])
  
}

# Explore the different components of the list
str(counts_ls_broad)



## Initiate empty list
counts_ls_specific <- list()

for (i in 1:length(cluster_names_specific)) {
  
  ## Extract indexes of columns in the global matrix that match a given cluster
  column_idx <- which(tstrsplit(colnames(aggr_counts_specific), "_")[[1]] == cluster_names_specific[i])
  
  ## Store corresponding sub-matrix as one element of a list
  counts_ls_specific[[i]] <- aggr_counts_specific[, column_idx]
  names(counts_ls_specific)[i] <- as.character(cluster_names_specific[i])
  
}

# Explore the different components of the list
str(counts_ls_specific)


# Reminder: explore structure of metadata
head(colData(sce))

# Extract sample-level variables
metadata <- colData(sce) %>% 
  as.data.frame() %>% 
  dplyr::select(TimeSurvival,Patient,sample_id)

dim(metadata)
head(metadata)

# Exclude duplicated rows
metadata <- metadata[!duplicated(metadata), ]

dim(metadata)
head(metadata)

# Rename rows
rownames(metadata) <- metadata$sample_id
head(metadata)


# Number of cells per sample and cluster
t <- table(colData(sce)$sample_id,colData(sce)$manualAnno_specific)
t2 <- table(colData(sce)$sample_id,colData(sce)$manualAnno_specific2)

t ### CANT COMPARE due to low count: Megakaryocte, Neutrophil, stromal
B cell CD34  CD4  CD8   DC Differentiating Stem Cell  EPC  gdT MAIT Megakaryocyte Monocyte Neutrophil   NK Plasmablast Stromal
113446_Short_Post    174  818  325  105  270                        48   98   64    4             0       75          8   78           7       0
113446_Short_Pre     204   30 1190  215  286                        12   22   97   16             0       18          0  123           7       0
114152_Long_Post      13  235  611  581    2                        21    8  246   67             0        8          0  113          42       0
114152_Long_Pre        1   12   13   26    1                         1    3   11    5             0        0          0    2           3       0
115585_Long_Post      24   40  143   32   13                        34  333    8    1             7       30          1   23          17       1
115585_Long_Pre      105   96  307  108   31                        59   38   43   13             2      522          8  125          14       0
118877_Long_Post      55   85  790  698   54                        35   12  236   24             5      137          2  403           8       2
118877_Long_Pre      401  345  923  277   13                       270   44   58   12             1       21          1  424          88       4
119325_Long_Post      51  313  444   85   11                        19  117   26    1             0       32         17   69          23       3
119325_Long_Pre       42  705  223  121   27                        31  128   59   14             1      195          6   44          15       0
119837_Long_Post      15   65  187  158    2                        29  491   98   28             4       41         16   21          11       2
119837_Long_Pre       10 1032  106  146   50                       114  190   77  112            19     1337        165   38           9       1
121015_Long_Post      15  106  833  707   28                        46   25  110    8             8      165        180  171           1       0
121015_Long_Pre       22  397  631  221   29                        50   24   33    7             4      140        115   47           4       0
121474_Long_Post      69  555  997  432    1                       336  494   61   84             4        4          2  316         157       0
121474_Long_Pre       30 1622  264   78   47                       196  521   15   23             2       69          7  439          15       0
121535_Long_Post      57  513  503  746    3                        48  855  161   69             3        7          1  210          56       4
121535_Long_Pre      129 1349  149  138   22                        71 1269   36   22             0       76          0   38         261      12
122739_Long_Post     283  538  504  532  162                       145   90  217   21             3      972          3  471         102       2
122739_Long_Pre      370  428 2070  641  147                       175   55  247   14             7      677          3  421          72       2
122759_Short_Post     53  233 1548  199    0                        23    0   19    0             0       44          0  140           6       0
122759_Short_Pre     124   82 1705  212    3                        16    0   17    2             2       87          7  193          14       0
122976_Long_Post       9   72 1627  507   22                        21   25  150    1             0       82          0  402           1       0
122976_Long_Pre       58   48  719  181   18                         9   23   55    0             4       90          0  161           1       0
24780_Short_Post      27  356  173   61    4                        58   32   49    8             1        5          1   30           7       1
24780_Short_Pre      187  559 1212  489   11                       256   88  149   60             0       19          1  145          13       0

t2 ## NOTE THAT IF YOU DO CELL # CUTOFF OF 10, ONLY FEW CLUSTERS HAVE SAMPLES THAT WORK FOR DEGS: CD4 NAIVE, CD4 TCM, CD56dim NK, CD8 TCM, CD8 TEM, CD8 TEMRA, Classical Mono, EPC, GMP, HSC, LMMP, mDC, Memory B, MEP, Vd1 gdt 
CD4 CTL CD4 Naive CD4 TCM CD4 TEM CD4 Treg CD56bright NK CD56dim NK CD8 Naive CD8 TCM CD8 TEM CD8 TEMRA Classical Mono Differentiating Stem Cell  EPC  GMP  HSC Intermediate Mono LMPP MAIT  mDC Megakaryocyte Memory B  MEP Naive B Neutrophil Non-classical Mono  pDC Plasmablast Stromal Vd1 gdT Vd2 gdT
113446_Short_Post       1       123     175       0       26             9         69         9      12      20        64             67                        48   98    1   10                 8  804    4  269             0      156    3      18          8                  0    1           7       0      59       5
113446_Short_Pre       15       575     550       2       48             8        115        14      22      53       126             13                        12   22    0    1                 3   29   16  281             0      172    0      32          0                  2    5           7       0      81      16
114152_Long_Post       43       113     339       0      116            10        103        41      18     234       288              6                        21    8    4  134                 0   71   67    0             0       12   28       1          0                  0    2          42       0     221      25
114152_Long_Pre         0         6       6       0        1             1          1         1       1      16         8              0                         1    3    0   11                 0    1    5    0             0        0    0       1          0                  0    1           3       0      11       0
115585_Long_Post        0        71      69       1        2            16          7         1       5      11        15             22                        34  333    7   16                 1    9    1   11             7       21   15       3          1                  0    2          17       1       8       0
115585_Long_Pre         9       143     150       1        4            39         86         3      15      42        48            419                        59   38  107   41                 4   36   13   27             2       97    9       8          8                  2    4          14       0      42       1
118877_Long_Post       93       404     256       1       36            10        393        13      24     263       398             85                        35   12   38   27                 6   48   24   22             5       45    9      10          2                  9   32           8       2     235       1
118877_Long_Pre        54       457     407       1        4            26        398        37      11      87       142             11                       270   44   22   50                 1  266   12    6             1      348   14      53          1                  2    7          88       4      50       8
119325_Long_Post        5       308     114      10        7             1         68         4       4       9        68             19                        19  117   22    5                 3  228    1    6             0       18   68      33         17                  0    5          23       3      26       0
119325_Long_Pre         3       113      98       3        6             4         40         1      12      44        64             68                        31  128  239    6                15  411   14   15             1       18  161      24          6                  0   12          15       0      57       2
119837_Long_Post       29        70      76       3        9             4         17         4       7      32       115             19                        29  491    3   10                11   24   28    2             4       12   30       3         16                  9    0          11       2      10      88
119837_Long_Pre         6        20      75       2        3            25         13         3      11      81        51            603                       114  190  871  144                 6  646  112   43            19       10   87       0        165                 12    7           9       1      15      62
121015_Long_Post       15       299     499       0       20            17        154        14      19     218       456            126                        46   25   38   26                15   43    8   21             8        8   15       7        180                  8    7           1       0     103       7
121015_Long_Pre         3       230     380       0       18             9         38         6      13      59       143             84                        50   24  128  142                 1  147    7   22             4       20   34       2        115                  1    7           4       0      32       1
121474_Long_Post        6       299     602       5       85            97        219        38      28     277        89              2                       336  494   15   99                 1  120   84    0             4       54  322      15          2                  0    1         157       0      19      42
121474_Long_Pre         0        69     181       2       12           242        197         2       4      54        18             34                       196  521  121  744                 9  510   23   17             2       21  264       9          7                  9   30          15       0       4      11
121535_Long_Post       74        34     330      54       11            36        174        26      94     336       290              3                        48  855   58    9                 1  293   69    2             3       47  153      10          1                  3    1          56       4     122      39
121535_Long_Pre        10        17     122       0        0            10         28         7      40      69        22             23                        71 1269  483   27                 1  797   22   21             0      121   68       8          0                 26    1         261      12      27       9
122739_Long_Post       24        47     397      13       23            62        409         5      22     150       355            787                       145   90   45  272               132  182   21   97             3      200   67      83          3                 25   65         102       2     202      15
122739_Long_Pre        33       791    1211       9       26            86        335        23      24     146       448            541                       175   55   60  147                76  205   14   98             7      266   34     104          3                 42   49          72       2     230      17
122759_Short_Post       9       709     710      24       96            13        127         6      67      51        75              1                        23    0    1   98                 3   40    0    0             0       43   94      10          0                 40    0           6       0      18       1
122759_Short_Pre        3      1131     532       8       31            22        171        18      99      63        32             17                        16    0    0   38                 4   18    2    3             2       47   26      77          7                 66    0          14       0      15       2
122976_Long_Post        5       452    1136       2       32            45        357       132     115      98       162             67                        21   25    9    9                 6   32    1    7             0        7   29       2          0                  2   15           1       0      93      57
122976_Long_Pre         6       251     448       1       13             8        153        24      37      22        98             68                         9   23    5   18                 4   13    0    8             4       44   13      14          0                 17   10           1       0      39      16
24780_Short_Post       30        21     117       2        3             0         30         2       2      17        40              1                        58   32   73  149                 0  119    8    2             1       25   19       2          1                  0    2           7       1      18      31
24780_Short_Pre        74       263     857       4       14             6        139        32      22     206       229              1                       256   88  151  109                 0  229   60    7             0      180   88       7          1                  0    4          13       0      56      93
###NOTE THAT LOT OF CELLS ARE LOW ABUNDANCE PER PATIENT PER TIMEPOINT


# Creating metadata list

## Initiate empty list
metadata_ls_broad <- list()

for (i in c("CD34")) {
  
  ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
  df <- data.frame(cluster_sample_id = colnames(counts_ls_broad[[i]]))
  
  ## Use tstrsplit() to separate cluster (cell type) and sample IDs
  df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
  df$sample_id  <- paste(tstrsplit(df$cluster_sample_id, "_")[[2]],tstrsplit(df$cluster_sample_id, "_")[[3]],tstrsplit(df$cluster_sample_id, "_")[[4]],sep="_")
  
  
  ## Retrieve cell count information for this cluster from global cell count table
  idx <- which(colnames(t) == unique(df$cluster_id))
  cell_counts <- t[, idx]
  
  ## Remove samples with less than 10 cell contributing to the cluster
  cell_counts <- cell_counts[cell_counts > 10]
  
  ## Match order of cell_counts and sample_ids
  sample_order <- match(df$sample_id, names(cell_counts))
  cell_counts <- cell_counts[sample_order]
  
  ## Append cell_counts to data frame
  df$cell_count <- cell_counts
  
  
  ## Join data frame (capturing metadata specific to cluster) to generic metadata
  df <- plyr::join(df, metadata, 
                   by = intersect(names(df), names(metadata)))
  
  ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
  rownames(df) <- df$cluster_sample_id
  
  #remove samples that don't pass cell # cutoff
  df <- df[complete.cases(df),]
  
  ## Store complete metadata for cluster i in list
  metadata_ls_broad[[i]] <- df
}

# Explore the different components of the list
str(metadata_ls_broad)




## Initiate empty list
metadata_ls_specific <- list()

for (i in 1:length(counts_ls_specific)) {
  
  ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
  df <- data.frame(cluster_sample_id = colnames(counts_ls_specific[[i]]))
  
  ## Use tstrsplit() to separate cluster (cell type) and sample IDs
  df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
  df$sample_id  <- paste(tstrsplit(df$cluster_sample_id, "_")[[2]],tstrsplit(df$cluster_sample_id, "_")[[3]],tstrsplit(df$cluster_sample_id, "_")[[4]],sep="_")
  
  
  ## Retrieve cell count information for this cluster from global cell count table
  idx <- which(colnames(t2) == unique(df$cluster_id))
  cell_counts <- t2[, idx]
  
  ## Remove samples with < 10 cell contributing to the cluster
  cell_counts <- cell_counts[cell_counts > 10]
  
  ## Match order of cell_counts and sample_ids
  sample_order <- match(df$sample_id, names(cell_counts))
  cell_counts <- cell_counts[sample_order]
  
  ## Append cell_counts to data frame
  df$cell_count <- cell_counts
  
  
  ## Join data frame (capturing metadata specific to cluster) to generic metadata
  df <- plyr::join(df, metadata, 
                   by = intersect(names(df), names(metadata)))
  
  ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
  rownames(df) <- df$cluster_sample_id

  ## Store complete metadata for cluster i in list
  metadata_ls_specific[[i]] <- df
  names(metadata_ls_specific)[i] <- unique(df$cluster_id)
  metadata_ls_specific[[i]] <- metadata_ls_specific[[i]][complete.cases(metadata_ls_specific[[i]]),]
  print(unique(df$cluster_id))
  print(table(metadata_ls_specific[[i]]$TimeSurvival))
}
[1] "GMP"

Long_Post   Long_Pre Short_Post  Short_Pre 
6          8          1          1 
[1] "CD4 TCM"

Long_Post   Long_Pre Short_Post  Short_Pre 
10          9          3          3 
[1] "Classical Mono"

Long_Post   Long_Pre Short_Post  Short_Pre 
7          9          1          2 
[1] "EPC"

Long_Post   Long_Pre Short_Post  Short_Pre 
9          9          2          2 
[1] "Differentiating Stem Cell"

Long_Post   Long_Pre Short_Post  Short_Pre 
10          8          3          3 
[1] "Plasmablast"

Long_Post  Long_Pre Short_Pre 
7         6         2 
[1] "Memory B"

Long_Post   Long_Pre Short_Post  Short_Pre 
8          8          3          3 
[1] "Neutrophil"

Long_Post  Long_Pre 
3         2 
[1] "Intermediate Mono"

Long_Post  Long_Pre 
3         2 
[1] "mDC"

Long_Post   Long_Pre Short_Post  Short_Pre 
4          7          1          1 
[1] "Non-classical Mono"

Long_Post   Long_Pre Short_Post  Short_Pre 
1          4          1          1 
[1] "CD8 TEMRA"

Long_Post   Long_Pre Short_Post  Short_Pre 
10          9          3          3 
[1] "LMPP"

Long_Post   Long_Pre Short_Post  Short_Pre 
9          9          3          3 
[1] "HSC"

Long_Post   Long_Pre Short_Post  Short_Pre 
6          9          2          2 
[1] "MEP"

Long_Post   Long_Pre Short_Post  Short_Pre 
9          8          2          2 
[1] "CD4 Naive"

Long_Post   Long_Pre Short_Post  Short_Pre 
10          9          3          3 
[1] "CD56dim NK"

Long_Post   Long_Pre Short_Post  Short_Pre 
9          9          3          3 
[1] "Naive B"

Long_Post   Long_Pre Short_Post  Short_Pre 
3          4          1          2 
[1] "Megakaryocyte"

Long_Pre 
1 
[1] "Vd1 gdT"

Long_Post   Long_Pre Short_Post  Short_Pre 
8          9          3          3 
[1] "CD8 TEM"

Long_Post   Long_Pre Short_Post  Short_Pre 
9         10          3          3 
[1] "pDC"

Long_Post  Long_Pre 
3         3 
[1] "CD56bright NK"

Long_Post   Long_Pre Short_Post  Short_Pre 
6          5          1          1 
[1] "CD4 Treg"

Long_Post   Long_Pre Short_Post  Short_Pre 
7          4          2          3 
[1] "CD4 CTL"

Long_Post   Long_Pre Short_Post  Short_Pre 
6          2          1          2 
[1] "CD8 Naive"

Long_Post  Long_Pre Short_Pre 
6         3         3 
[1] "CD8 TCM"

Long_Post   Long_Pre Short_Post  Short_Pre 
7          8          2          3 
[1] "Vd2 gdT"

Long_Post   Long_Pre Short_Post  Short_Pre 
6          4          1          2 
[1] "MAIT"

Long_Post  Long_Pre Short_Pre 
6         7         2 
[1] "Stromal"

Long_Pre 
1 
[1] "CD4 TEM"

Long_Post Short_Post 
2          1 




# Explore the different components of the list
str(metadata_ls_specific)



######## Creating a DESeq2 object #############

DESeq2_results_broad <-data.frame()
DEGs_compiled_broad <- data.frame()

for (broadcelltype in c("CD34")){
  cluster_counts <- counts_ls_broad[[broadcelltype]]
  cluster_metadata <- metadata_ls_broad[[broadcelltype]]
  cluster_metadata$Patient <- factor(cluster_metadata$Patient)
  
  cluster_counts <- cluster_counts[,cluster_metadata$cluster_sample_id]
  
  all(colnames(cluster_counts) == rownames(cluster_metadata))
  head(cluster_metadata)
  
  # Create DESeq2 object    
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ TimeSurvival)
  keep <- rowSums(counts(dds)) > 10 #filter non expressed genes 
  dds <- dds[keep,]#filter lowly expressed genes
  dds <- DESeq(dds)
  
  # Transform counts for data visualization
  rld <- rlog(dds, blind=TRUE)
  
  # Plot PCA
  DESeq2::plotPCA(rld, ntop = 500, intgroup = "TimeSurvival") +ggtitle(broadcelltype)  
  
  # Extract the rlog matrix from the object and compute pairwise correlation values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  
  # Plot heatmap
  pheatmap(rld_cor, annotation = cluster_metadata[, c("Patient","TimeSurvival"), drop=F])
  
  
  # Generate results object
  res_short_pre_post <- results(dds,contrast=c("TimeSurvival","Short_Pre","Short_Post"),alpha=0.05)
  res_long_pre_post <- results(dds,contrast=c("TimeSurvival","Long_Pre","Long_Post"),alpha=0.05)
  res_short_long_pre <- results(dds,contrast=c("TimeSurvival","Short_Pre","Long_Pre"),alpha=0.05)
  res_short_long_post <- results(dds,contrast=c("TimeSurvival","Short_Post","Long_Post"),alpha=0.05)
  
  
  # Turn the DESeq2 results object into a tibble for use with tidyverse functions
  res_tbl_short_pre_post  <- res_short_pre_post  %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    arrange(stat)
  
  sig_res_short_pre_post <- dplyr::filter(res_tbl_short_pre_post, padj < 0.05) %>%
    dplyr::arrange(padj)
  
  
  
  res_tbl_long_pre_post  <- res_long_pre_post  %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    arrange(stat)
  
  sig_res_long_pre_post <- dplyr::filter(res_tbl_long_pre_post, padj < 0.05) %>%
    dplyr::arrange(padj)
  
  
  
  res_tbl_short_long_pre  <- res_short_long_pre  %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    arrange(stat)
  
  sig_res_short_long_pre <- dplyr::filter(res_tbl_short_long_pre, padj < 0.05) %>%
    dplyr::arrange(padj)
  
  
  res_tbl_short_long_post  <- res_short_long_post  %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    arrange(stat)
  
  sig_res_short_long_post <- dplyr::filter(res_tbl_short_long_post, padj < 0.05) %>%
    dplyr::arrange(padj)
  
  
  res_tbl_short_pre_post$Celltype <- broadcelltype
  res_tbl_short_pre_post$Comparison <- "ShortPre_ShortPost"
  res_tbl_long_pre_post$Celltype <- broadcelltype
  res_tbl_long_pre_post$Comparison <- "LongPre_LongPost"
  res_tbl_short_long_pre$Celltype <- broadcelltype
  res_tbl_short_long_pre$Comparison <- "ShortPre_LongPre"
  res_tbl_short_long_post$Celltype <- broadcelltype
  res_tbl_short_long_post$Comparison <- "ShortPost_LongPost"
  
  sig_res_short_pre_post$Celltype <- broadcelltype
  sig_res_short_pre_post$Comparison <- "ShortPre_ShortPost"
  sig_res_long_pre_post$Celltype <- broadcelltype
  sig_res_long_pre_post$Comparison <- "LongPre_LongPost"
  sig_res_short_long_pre$Celltype <- broadcelltype
  sig_res_short_long_pre$Comparison <- "ShortPre_LongPre"
  sig_res_short_long_post$Celltype <- broadcelltype
  sig_res_short_long_post$Comparison <- "ShortPost_LongPost"
  
  DESeq2_results_broad <- rbind(DESeq2_results_broad,res_tbl_short_pre_post,res_tbl_long_pre_post,res_tbl_short_long_pre,res_tbl_short_long_post)
  DEGs_compiled_broad <- rbind(DEGs_compiled_broad,sig_res_short_pre_post,sig_res_long_pre_post,sig_res_short_long_pre,sig_res_short_long_post)
  
  print(broadcelltype)
  print("sig_res_short_pre_post")
  print(sig_res_short_pre_post)
  print("sig_res_long_pre_post")
  print(sig_res_long_pre_post)
  print("sig_res_short_long_pre")
  print(sig_res_short_long_pre)
  print("sig_res_short_long_post")
  print(sig_res_short_long_post)
}


DESeq2_results_specific <-data.frame()
DEGs_compiled_specific <- data.frame()

for (specificcelltype in c("CD4 TCM","EPC","Memory B","CD8 TEMRA","LMPP","HSC","MEP","CD4 Naive", "CD56dim NK","Vd1 gdT","CD8 TEM", "CD4 Treg","CD8 TCM")){
  tryCatch({ print(specificcelltype)
    cluster_counts <- counts_ls_specific[[specificcelltype]]
    cluster_metadata <- metadata_ls_specific[[specificcelltype]]
    cluster_metadata$Patient <- factor(cluster_metadata$Patient)
    
    cluster_counts <- cluster_counts[,cluster_metadata$cluster_sample_id]
    
    all(colnames(cluster_counts) == rownames(cluster_metadata))
    head(cluster_metadata)
    
    # Create DESeq2 object    
    dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                  colData = cluster_metadata, 
                                  design = ~ TimeSurvival)
    keep <- rowSums(counts(dds)) >= 10 #filter non expressed genes 
    dds <- dds[keep,]#filter lowly expressed genes
    dds <- DESeq(dds)
    
    # Transform counts for data visualization
    rld <- rlog(dds, blind=TRUE)
    
    # Plot PCA
    DESeq2::plotPCA(rld, ntop = 500, intgroup = "TimeSurvival") +ggtitle(specificcelltype)  
    
    # Extract the rlog matrix from the object and compute pairwise correlation values
    rld_mat <- assay(rld)
    rld_cor <- cor(rld_mat)
    
    # Plot heatmap
    pheatmap(rld_cor, annotation = cluster_metadata[, c("Patient","TimeSurvival"), drop=F])
    
    
    # Generate results object
    res_short_pre_post <- results(dds,contrast=c("TimeSurvival","Short_Pre","Short_Post"),alpha=0.05)
    res_long_pre_post <- results(dds,contrast=c("TimeSurvival","Long_Pre","Long_Post"),alpha=0.05)
    res_short_long_pre <- results(dds,contrast=c("TimeSurvival","Short_Pre","Long_Pre"),alpha=0.05)
    res_short_long_post <- results(dds,contrast=c("TimeSurvival","Short_Post","Long_Post"),alpha=0.05)
    
    
    # Turn the DESeq2 results object into a tibble for use with tidyverse functions
    res_tbl_short_pre_post  <- res_short_pre_post  %>%
      data.frame() %>%
      rownames_to_column(var = "gene") %>%
      as_tibble() %>%
      arrange(stat)
    
    sig_res_short_pre_post <- dplyr::filter(res_tbl_short_pre_post, padj < 0.05) %>%
      dplyr::arrange(padj)
    
    
    
    res_tbl_long_pre_post  <- res_long_pre_post  %>%
      data.frame() %>%
      rownames_to_column(var = "gene") %>%
      as_tibble() %>%
      arrange(stat)
    
    sig_res_long_pre_post <- dplyr::filter(res_tbl_long_pre_post, padj < 0.05) %>%
      dplyr::arrange(padj)
    
    
    
    res_tbl_short_long_pre  <- res_short_long_pre  %>%
      data.frame() %>%
      rownames_to_column(var = "gene") %>%
      as_tibble() %>%
      arrange(stat)
    
    sig_res_short_long_pre <- dplyr::filter(res_tbl_short_long_pre, padj < 0.05) %>%
      dplyr::arrange(padj)
    
    
    res_tbl_short_long_post  <- res_short_long_post  %>%
      data.frame() %>%
      rownames_to_column(var = "gene") %>%
      as_tibble() %>%
      arrange(stat)
    
    sig_res_short_long_post <- dplyr::filter(res_tbl_short_long_post, padj < 0.05) %>%
      dplyr::arrange(padj)
    
    
    res_tbl_short_pre_post$Celltype <- specificcelltype
    res_tbl_short_pre_post$Comparison <- "ShortPre_ShortPost"
    res_tbl_long_pre_post$Celltype <- specificcelltype
    res_tbl_long_pre_post$Comparison <- "LongPre_LongPost"
    res_tbl_short_long_pre$Celltype <- specificcelltype
    res_tbl_short_long_pre$Comparison <- "ShortPre_LongPre"
    res_tbl_short_long_post$Celltype <- specificcelltype
    res_tbl_short_long_post$Comparison <- "ShortPost_LongPost"
    
    sig_res_short_pre_post$Celltype <- specificcelltype
    sig_res_short_pre_post$Comparison <- "ShortPre_ShortPost"
    sig_res_long_pre_post$Celltype <- specificcelltype
    sig_res_long_pre_post$Comparison <- "LongPre_LongPost"
    sig_res_short_long_pre$Celltype <- specificcelltype
    sig_res_short_long_pre$Comparison <- "ShortPre_LongPre"
    sig_res_short_long_post$Celltype <- specificcelltype
    sig_res_short_long_post$Comparison <- "ShortPost_LongPost"
    
    DESeq2_results_specific <- rbind(DESeq2_results_specific,res_tbl_short_pre_post,res_tbl_long_pre_post,res_tbl_short_long_pre,res_tbl_short_long_post)
    DEGs_compiled_specific <- rbind(DEGs_compiled_specific,sig_res_short_pre_post,sig_res_long_pre_post,sig_res_short_long_pre,sig_res_short_long_post)
    
    print(specificcelltype)
    print("sig_res_short_pre_post")
    print(sig_res_short_pre_post)
    print("sig_res_long_pre_post")
    print(sig_res_long_pre_post)
    print("sig_res_short_long_pre")
    print(sig_res_short_long_pre)
    print("sig_res_short_long_post")
    print(sig_res_short_long_post)
  }, 
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})}


for (specificcelltype in c("GMP","Neutrophil","Intermediate Mono","mDC","pDC","CD56bright NK")){ # Not enough samples for short pre or short post so need to exclude it from code
  tryCatch({print(specificcelltype)
    cluster_counts <- counts_ls_specific[[specificcelltype]]
    cluster_metadata <- metadata_ls_specific[[specificcelltype]]
    cluster_metadata$Patient <- factor(cluster_metadata$Patient)
    
    cluster_counts <- cluster_counts[,cluster_metadata$cluster_sample_id]
    
    all(colnames(cluster_counts) == rownames(cluster_metadata))
    head(cluster_metadata)
    
    # Create DESeq2 object    
    dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                  colData = cluster_metadata, 
                                  design = ~ TimeSurvival)
    keep <- rowSums(counts(dds)) >= 10 #filter non expressed genes 
    dds <- dds[keep,]#filter lowly expressed genes
    dds <- DESeq(dds)
    
    # Transform counts for data visualization
    rld <- rlog(dds, blind=TRUE)
    
    # Plot PCA
    DESeq2::plotPCA(rld, ntop = 500, intgroup = "TimeSurvival") +ggtitle(specificcelltype)  
    
    # Extract the rlog matrix from the object and compute pairwise correlation values
    rld_mat <- assay(rld)
    rld_cor <- cor(rld_mat)
    
    # Plot heatmap
    pheatmap(rld_cor, annotation = cluster_metadata[, c("Patient","TimeSurvival"), drop=F])
    
    # Generate results object
    res_long_pre_post <- results(dds,contrast=c("TimeSurvival","Long_Pre","Long_Post"),alpha=0.05)
    
    res_tbl_long_pre_post  <- res_long_pre_post  %>%
      data.frame() %>%
      rownames_to_column(var = "gene") %>%
      as_tibble() %>%
      arrange(stat)
    
    sig_res_long_pre_post <- dplyr::filter(res_tbl_long_pre_post, padj < 0.05) %>%
      dplyr::arrange(padj)
    
    res_tbl_long_pre_post$Celltype <- specificcelltype
    res_tbl_long_pre_post$Comparison <- "LongPre_LongPost"
   
    sig_res_long_pre_post$Celltype <- specificcelltype
    sig_res_long_pre_post$Comparison <- "LongPre_LongPost"
    
    DESeq2_results_specific <- rbind(DESeq2_results_specific,res_tbl_long_pre_post)
    DEGs_compiled_specific <- rbind(DEGs_compiled_specific,sig_res_long_pre_post)
    
    print(specificcelltype)
   
    print("sig_res_long_pre_post")
    print(sig_res_long_pre_post)
  }, 
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})}



for (specificcelltype in c("Classical Mono","Plasmablast","Naive B","CD4 CTL", "CD8 Naive","Vd2 gdT","MAIT")){ # Not enough samples for short post so need to exclude it from code
  tryCatch({print(specificcelltype)
    cluster_counts <- counts_ls_specific[[specificcelltype]]
    cluster_metadata <- metadata_ls_specific[[specificcelltype]]
    cluster_metadata$Patient <- factor(cluster_metadata$Patient)
    
    cluster_counts <- cluster_counts[,cluster_metadata$cluster_sample_id]
    
    all(colnames(cluster_counts) == rownames(cluster_metadata))
    head(cluster_metadata)
    
    # Create DESeq2 object    
    dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                  colData = cluster_metadata, 
                                  design = ~ TimeSurvival)
    keep <- rowSums(counts(dds)) >= 10 #filter non expressed genes 
    dds <- dds[keep,]#filter lowly expressed genes
    dds <- DESeq(dds)
    
    # Transform counts for data visualization
    rld <- rlog(dds, blind=TRUE)
    
    # Plot PCA
    DESeq2::plotPCA(rld, ntop = 500, intgroup = "TimeSurvival") +ggtitle(specificcelltype)  
    
    # Extract the rlog matrix from the object and compute pairwise correlation values
    rld_mat <- assay(rld)
    rld_cor <- cor(rld_mat)
    
    # Plot heatmap
    pheatmap(rld_cor, annotation = cluster_metadata[, c("Patient","TimeSurvival"), drop=F])
    
    
    # Generate results object
    res_long_pre_post <- results(dds,contrast=c("TimeSurvival","Long_Pre","Long_Post"),alpha=0.05)
    res_short_long_pre <- results(dds,contrast=c("TimeSurvival","Short_Pre","Long_Pre"),alpha=0.05)

    res_tbl_long_pre_post  <- res_long_pre_post  %>%
      data.frame() %>%
      rownames_to_column(var = "gene") %>%
      as_tibble() %>%
      arrange(stat)
    
    sig_res_long_pre_post <- dplyr::filter(res_tbl_long_pre_post, padj < 0.05) %>%
      dplyr::arrange(padj)
    
     res_tbl_short_long_pre  <- res_short_long_pre  %>%
    data.frame() %>%
     rownames_to_column(var = "gene") %>%
     as_tibble() %>%
     arrange(stat)
    
     sig_res_short_long_pre <- dplyr::filter(res_tbl_short_long_pre, padj < 0.05) %>%
    dplyr::arrange(padj)
    

    res_tbl_long_pre_post$Celltype <- specificcelltype
    res_tbl_long_pre_post$Comparison <- "LongPre_LongPost"
    res_tbl_short_long_pre$Celltype <- specificcelltype
    res_tbl_short_long_pre$Comparison <- "ShortPre_LongPre"

    sig_res_long_pre_post$Celltype <- specificcelltype
    sig_res_long_pre_post$Comparison <- "LongPre_LongPost"
    sig_res_short_long_pre$Celltype <- specificcelltype
    sig_res_short_long_pre$Comparison <- "ShortPre_LongPre"

    DESeq2_results_specific <- rbind(DESeq2_results_specific,res_tbl_long_pre_post)
    DEGs_compiled_specific <- rbind(DEGs_compiled_specific,sig_res_long_pre_post)
    
    print(specificcelltype)
    print("sig_res_long_pre_post")
    print(sig_res_long_pre_post)
     print("sig_res_short_long_pre")
     print(sig_res_short_long_pre)
  }, 
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})}


#### Perform GSEA analysis on the DEGs (https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/scRNAseq_workshop_3.html, https://github.com/ctlab/fgsea/issues/50) ######

hallmark<- msigdbr(species = "Homo sapiens", category = "H")

bpGO <- msigdbr(species = "Homo sapiens", category = "C5",subcategory = "GO:BP")
mfGO <- msigdbr(species = "Homo sapiens", category = "C5",subcategory = "GO:MF")
ccGO <- msigdbr(species = "Homo sapiens", category = "C5",subcategory = "GO:CC")


fgsea_sets_hallmark<- hallmark %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets_bpGO<- bpGO %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets_mfGO<- mfGO %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets_ccGO<- ccGO %>% split(x = .$gene_symbol, f = .$gs_name)

fgsea_SENMAYO <- list(c("ACVR1B","ANG","ANGPT1","ANGPTL4","AREG","AXL","BEX3","BMP2","BMP6","C3","CCL1","CCL13","CCL16","CCL2","CCL20","CCL24","CCL26","CCL3","CCL3L1","CCL4","CCL5","CCL7","CCL8","CD55","CD9","CSF1","CSF2","CSF2RB","CST4","CTNNB1","CTSB","CXCL1","CXCL10","CXCL12","CXCL16","CXCL2","CXCL3","CXCL8","CXCR2","DKK1","EDN1","EGF","EGFR","EREG","ESM1","ETS2","FAS","FGF1","FGF2","FGF7","GDF15","GEM","GMFG","HGF","HMGB1","ICAM1","ICAM3","IGF1","IGFBP1","IGFBP2","IGFBP3","IGFBP4","IGFBP5","IGFBP6","IGFBP7","IL10","IL13","IL15","IL18","IL1A","IL1B","IL2","IL32","IL6","IL6ST","IL7","INHA","IQGAP2","ITGA2","ITPKA","JUN","KITLG","LCP1","MIF","MMP1","MMP10","MMP12","MMP13","MMP14","MMP2","MMP3","MMP9","NAP1L4","NRG1","PAPPA","PECAM1","PGF","PIGF","PLAT","PLAU","PLAUR","PTBP1","PTGER2","PTGES","RPS6KA5","SCAMP4","SELPLG","SEMA3F","SERPINB4","SERPINE1","SERPINE2","SPP1","SPX","TIMP2","TNF","TNFRSF10C","TNFRSF11B","TNFRSF1A","TNFRSF1B","TUBGCP2","VEGFA","VEGFC","VGF","WNT16","WNT2"))
names(fgsea_SENMAYO) <- "SenMayo"



df_GSEA_broad_DESeq2 <- data.frame()


#Hallmark

for (i in unique(DESeq2_results_broad$Comparison))
{for (x in unique(DESeq2_results_broad[DESeq2_results_broad$Comparison == i,]$Celltype))
{print(paste(i,x,sep = " | "))
  
  fgsea_ranks <- DESeq2_results_broad[DESeq2_results_broad$Comparison == i & DESeq2_results_broad$Celltype == x,c("gene","stat")]
  row.names(fgsea_ranks) <- fgsea_ranks$gene
  fgsea_ranks <-fgsea_ranks[order(-fgsea_ranks$stat),]
  fgsea_ranks <- deframe(fgsea_ranks)
  
  if (length(fgsea_ranks[fgsea_ranks != 0]) > 10) {
    fgseaRes <- fgseaMultilevel(pathways = fgsea_sets_hallmark,  stats = fgsea_ranks)
    
    fgseaResTidy <- fgseaRes %>%
      as_tibble() %>%
      arrange(desc(NES))
    
    fgseaResTidy %>% 
      dplyr::select(-leadingEdge, -ES) %>% 
      arrange(padj) %>% 
      DT::datatable()
    
    if (nrow(fgseaResTidy %>% filter(padj < 0.05)) > 0) {
      
      fgsea_df<-data.frame(fgseaResTidy %>% filter(padj < 0.05))
      fgsea_df$Comparison <- i
      fgsea_df$Celltype <- x
      fgsea_df$gsea_type <- "Hallmark"
      
      df_GSEA_broad_DESeq2 <- rbind(df_GSEA_broad_DESeq2, fgsea_df)} else {print("no significant Pathways")}
  } else {print("no DEGs")}}}



######## GSEA on Specific Anno comparison

df_GSEA_specific_DESeq2 <- data.frame()

#hallmark

for (i in unique(DESeq2_results_specific$Comparison))
{for (x in unique(DESeq2_results_specific[DESeq2_results_specific$Comparison == i,]$Celltype))
{print(paste(i,x,sep = " | "))
  
  fgsea_ranks <- DESeq2_results_specific[DESeq2_results_specific$Comparison == i & DESeq2_results_specific$Celltype == x,c("gene","stat")]
  fgsea_ranks <-fgsea_ranks[!duplicated(fgsea_ranks),]
  row.names(fgsea_ranks) <- fgsea_ranks$gene
  fgsea_ranks <-fgsea_ranks[order(-fgsea_ranks$stat),]
  fgsea_ranks <- deframe(fgsea_ranks)
  
  if (length(fgsea_ranks[fgsea_ranks != 0]) > 10) {
    fgseaRes <- fgseaMultilevel(pathways = fgsea_sets_hallmark,  stats = fgsea_ranks)
    
    fgseaResTidy <- fgseaRes %>%
      as_tibble() %>%
      arrange(desc(NES))
    
    fgseaResTidy %>% 
      dplyr::select(-leadingEdge, -ES) %>% 
      arrange(padj) %>% 
      DT::datatable()
    
    if (nrow(fgseaResTidy %>% filter(padj < 0.05)) > 0) {
      
      fgsea_df<-data.frame(fgseaResTidy %>% filter(padj < 0.05))
      fgsea_df$Comparison <- i
      fgsea_df$Celltype <- x
      fgsea_df$gsea_type <- "Hallmark"
      
      df_GSEA_specific_DESeq2 <- rbind(df_GSEA_specific_DESeq2, fgsea_df)} else {print("no significant Pathways")}
  } else {print("no DEGs")}}}


#SENMAYO

for (i in unique(DESeq2_results_specific$Comparison))
{for (x in unique(DESeq2_results_specific[DESeq2_results_specific$Comparison == i,]$Celltype))
{print(paste(i,x,sep = " | "))
  
  fgsea_ranks <- DESeq2_results_specific[DESeq2_results_specific$Comparison == i & DESeq2_results_specific$Celltype == x,c("gene","stat")]
  fgsea_ranks <-fgsea_ranks[!duplicated(fgsea_ranks),]
  row.names(fgsea_ranks) <- fgsea_ranks$gene
  fgsea_ranks <-fgsea_ranks[order(-fgsea_ranks$stat),]
  fgsea_ranks <- deframe(fgsea_ranks)
  
  if (length(fgsea_ranks[fgsea_ranks != 0]) > 10) {
    fgseaRes <- fgseaMultilevel(pathways = fgsea_SENMAYO,  stats = fgsea_ranks)
    
    fgseaResTidy <- fgseaRes %>%
      as_tibble() %>%
      arrange(desc(NES))
    
    fgseaResTidy %>% 
      dplyr::select(-leadingEdge, -ES) %>% 
      arrange(padj) %>% 
      DT::datatable()
    
    if (nrow(fgseaResTidy %>% filter(padj < 0.01)) > 0) {
      
      fgsea_df<-data.frame(fgseaResTidy %>% filter(padj < 0.01))
      fgsea_df$Comparison <- i
      fgsea_df$Celltype <- x
      fgsea_df$gsea_type <- "SenMayo"
      
      df_GSEA_specific_DESeq2 <- rbind(df_GSEA_specific_DESeq2, fgsea_df)} else {print("no significant Pathways")}
  } else {print("no DEGs")}}}






#PLot GSEA enrichment Senescence
fgsea_ranks <- DESeq2_results_specific[DESeq2_results_specific$Comparison == "ShortPre_LongPre" & DESeq2_results_specific$Celltype == "LMPP",c("gene","stat")]
fgsea_ranks <-fgsea_ranks[!duplicated(fgsea_ranks),]
row.names(fgsea_ranks) <- fgsea_ranks$gene
fgsea_ranks <-fgsea_ranks[order(-fgsea_ranks$stat),]
fgsea_ranks <- deframe(fgsea_ranks)
fgseaRes <- fgseaMultilevel(pathways = fgsea_SENMAYO,  stats = fgsea_ranks)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

plotEnrichment(fgsea_SENMAYO[["SenMayo"]],stats=fgsea_ranks)




####################################Make figures of DESeq2 pseudobulk GSEA RESULTS ###
##Pseudobulk CD34 

df.R<-data.frame(unique(c(df_GSEA_broad_DESeq2[df_GSEA_broad_DESeq2$gsea_type == "Hallmark",]$pathway)))
colnames(df.R) <-"pathway"
df.R<- merge(df.R,df_GSEA_broad_DESeq2[df_GSEA_broad_DESeq2$gsea_type == "Hallmark" & df_GSEA_broad_DESeq2$Comparison == "ShortPre_LongPre",c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,df_GSEA_broad_DESeq2[df_GSEA_broad_DESeq2$gsea_type == "Hallmark" & df_GSEA_broad_DESeq2$Comparison == "ShortPost_LongPost",c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,df_GSEA_broad_DESeq2[df_GSEA_broad_DESeq2$gsea_type == "Hallmark" & df_GSEA_broad_DESeq2$Comparison == "ShortPre_ShortPost",c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,df_GSEA_broad_DESeq2[df_GSEA_broad_DESeq2$gsea_type == "Hallmark" & df_GSEA_broad_DESeq2$Comparison == "LongPre_LongPost",c("pathway","NES")],by="pathway", all.x = T)
colnames(df.R)<-c("pathway","ShortPre_LongPre","ShortPost_LongPost","ShortPre_ShortPost","LongPre_LongPost")
df.R[is.na(df.R)]<-0
df.R$pathway_simple <- apply(df.R,1,function(x){unlist(strsplit(x[1],"MARK_"))[2]})


df.R2<-data.frame(unique(c(df_GSEA_broad_DESeq2[df_GSEA_broad_DESeq2$gsea_type == "Hallmark",]$pathway)))
colnames(df.R2) <-"pathway"
df.R2<- merge(df.R2,df_GSEA_broad_DESeq2[df_GSEA_broad_DESeq2$gsea_type == "Hallmark" & df_GSEA_broad_DESeq2$Comparison == "ShortPre_LongPre",c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,df_GSEA_broad_DESeq2[df_GSEA_broad_DESeq2$gsea_type == "Hallmark" & df_GSEA_broad_DESeq2$Comparison == "ShortPost_LongPost",c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,df_GSEA_broad_DESeq2[df_GSEA_broad_DESeq2$gsea_type == "Hallmark" & df_GSEA_broad_DESeq2$Comparison == "ShortPre_ShortPost",c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,df_GSEA_broad_DESeq2[df_GSEA_broad_DESeq2$gsea_type == "Hallmark" & df_GSEA_broad_DESeq2$Comparison == "LongPre_LongPost",c("pathway","padj")],by="pathway", all.x = T)
colnames(df.R2)<-c("pathway","ShortPre_LongPre","ShortPost_LongPost","ShortPre_ShortPost","LongPre_LongPost")
df.R2[is.na(df.R2)]<-1
df.R2$pathway_simple <- apply(df.R2,1,function(x){unlist(strsplit(x[1],"MARK_"))[2]})
df.R2[,c(2:5)]<- -log(df.R2[,c(2:5)])

df.R_M <- reshape2::melt(df.R[,c(6,2:5)])
df.R2_M <- reshape2::melt(df.R2[,c(6,2:5)])
df.R3_M <- merge(df.R_M,df.R2_M, by= c("pathway_simple","variable"))

df.R3_M <- df.R3_M[df.R3_M$pathway_simple %in% c("EPITHELIAL_MESENCHYMAL_TRANSITION", "MYOGENESIS","XENOBIOTIC_METABOLISM","HEME_METABOLISM","MTORC1_SIGNALING","TNFA_SIGNALING_VIA_NFKB","APOPTOSIS","INFLAMMATORY_RESPONSE","INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE","MYC_TARGETS_V1","MYC_TARGETS_V2","G2M_CHECKPOINT","E2F_TARGETS"),]

df.R3_M$variable <- factor(df.R3_M$variable, levels = c("ShortPre_LongPre","ShortPost_LongPost","ShortPre_ShortPost","LongPre_LongPost"))
df.R3_M$pathway_simple <- factor(df.R3_M$pathway_simple, levels = rev(c("EPITHELIAL_MESENCHYMAL_TRANSITION", "MYOGENESIS","XENOBIOTIC_METABOLISM","HEME_METABOLISM","MTORC1_SIGNALING","TNFA_SIGNALING_VIA_NFKB","APOPTOSIS","INFLAMMATORY_RESPONSE","INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE","MYC_TARGETS_V1","MYC_TARGETS_V2","G2M_CHECKPOINT","E2F_TARGETS"))) 
ggplot(df.R3_M, aes(x=variable, y = pathway_simple, color = value.x, size = value.y)) + 
  geom_point() + 
  scale_color_gradientn(colours =rev(c(brewer.pal(n = 5, name = "PuOr"))),limits=c(-4,4)) + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  ylab('') +
  theme(axis.ticks = element_blank()) 


######## SPECIFIC ANNO DESEQ2 PSEUDO GSEA ###

#Short Pretreatment vs Short Posttreatment
df.R<-data.frame(unique(c(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark" & df_GSEA_specific_DESeq2$Comparison == "ShortPre_ShortPost",]$pathway)))
colnames(df.R) <-"pathway"
for (i in unique(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark"& df_GSEA_specific_DESeq2$Comparison == "ShortPre_ShortPost",]$Celltype)) {
  print(i)
  df.R<- merge(df.R,df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark" & df_GSEA_specific_DESeq2$Comparison == "ShortPre_ShortPost" & df_GSEA_specific_DESeq2$Celltype == i,c("pathway","NES")],by="pathway", all.x = T)
}
colnames(df.R)<-c("pathway",c(unique(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark"& df_GSEA_specific_DESeq2$Comparison == "ShortPre_ShortPost",]$Celltype)))
df.R[is.na(df.R)]<-0
df.R$pathway_simple <- apply(df.R,1,function(x){unlist(strsplit(x[1],"MARK_"))[2]})


df.R2<-data.frame(unique(c(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark"& df_GSEA_specific_DESeq2$Comparison == "ShortPre_ShortPost",]$pathway)))
colnames(df.R2) <-"pathway"
for (i in unique(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark"& df_GSEA_specific_DESeq2$Comparison == "ShortPre_ShortPost",]$Celltype)) {
  print(i)
  df.R2<- merge(df.R2,df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark" & df_GSEA_specific_DESeq2$Comparison == "ShortPre_ShortPost" & df_GSEA_specific_DESeq2$Celltype == i,c("pathway","padj")],by="pathway", all.x = T)
}
colnames(df.R2)<-c("pathway",c(unique(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark"& df_GSEA_specific_DESeq2$Comparison == "ShortPre_ShortPost",]$Celltype)))
df.R2[is.na(df.R2)]<-1
df.R2$pathway_simple <- apply(df.R2,1,function(x){unlist(strsplit(x[1],"MARK_"))[2]})
df.R2[,c(2:13)]<- -log(df.R2[,c(2:13)])

df.R_M <- reshape2::melt(df.R[,c(14,2:13)])
df.R2_M <- reshape2::melt(df.R2[,c(14,2:13)])
df.R3_M <- merge(df.R_M,df.R2_M, by= c("pathway_simple","variable"))

df.R3_M$variable <- factor(df.R3_M$variable, levels = c("EPC","MEP","HSC","LMPP","Memory B","CD56dim NK", "Vd1 gdT","CD4 Naive","CD4 TCM","CD4 Treg","CD8 TCM","CD8 TEM","CD8 TEMRA"))
ggplot(df.R3_M, aes(x=variable, y = pathway_simple, color = value.x, size = value.y)) + 
  geom_point() + 
  scale_color_gradientn(colours =rev(c(brewer.pal(n = 5, name = "PuOr"))),limits=c(-4,4)) + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  ylab('') +
  theme(axis.ticks = element_blank()) 



#Long Pretreatment vs Long Posttreatment
df.R<-data.frame(unique(c(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark" & df_GSEA_specific_DESeq2$Comparison == "LongPre_LongPost",]$pathway)))
colnames(df.R) <-"pathway"
for (i in unique(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark"& df_GSEA_specific_DESeq2$Comparison == "LongPre_LongPost",]$Celltype)) {
  print(i)
  df.R<- merge(df.R,df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark" & df_GSEA_specific_DESeq2$Comparison == "LongPre_LongPost" & df_GSEA_specific_DESeq2$Celltype == i,c("pathway","NES")],by="pathway", all.x = T)
}
colnames(df.R)<-c("pathway",c(unique(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark"& df_GSEA_specific_DESeq2$Comparison == "LongPre_LongPost",]$Celltype)))
df.R[is.na(df.R)]<-0
df.R$pathway_simple <- apply(df.R,1,function(x){unlist(strsplit(x[1],"MARK_"))[2]})


df.R2<-data.frame(unique(c(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark"& df_GSEA_specific_DESeq2$Comparison == "LongPre_LongPost",]$pathway)))
colnames(df.R2) <-"pathway"
for (i in unique(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark"& df_GSEA_specific_DESeq2$Comparison == "LongPre_LongPost",]$Celltype)) {
  print(i)
  df.R2<- merge(df.R2,df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark" & df_GSEA_specific_DESeq2$Comparison == "LongPre_LongPost" & df_GSEA_specific_DESeq2$Celltype == i,c("pathway","padj")],by="pathway", all.x = T)
}
colnames(df.R2)<-c("pathway",c(unique(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark"& df_GSEA_specific_DESeq2$Comparison == "LongPre_LongPost",]$Celltype)))
df.R2[is.na(df.R2)]<-1
df.R2$pathway_simple <- apply(df.R2,1,function(x){unlist(strsplit(x[1],"MARK_"))[2]})
df.R2[,c(2:27)]<- -log(df.R2[,c(2:27)])

df.R_M <- reshape2::melt(df.R[,c(28,2:27)])
df.R2_M <- reshape2::melt(df.R2[,c(28,2:27)])
df.R3_M <- merge(df.R_M,df.R2_M, by= c("pathway_simple","variable"))

df.R3_M$variable <- factor(df.R3_M$variable, levels = c("EPC","MEP","HSC","LMPP","GMP","Classical Mono","Intermediate Mono","mDC","pDC","Naive B","Memory B","Plasmablast","CD56dim NK","CD56bright NK", "Vd1 gdT", "Vd2 gdT","CD4 Naive","CD4 TCM","CD4 CTL","CD4 Treg","CD8 Naive","CD8 TCM","CD8 TEM","CD8 TEMRA","MAIT","Neutrophil"))
ggplot(df.R3_M, aes(x=variable, y = pathway_simple, color = value.x, size = value.y)) + 
  geom_point() + 
  scale_color_gradientn(colours =rev(c(brewer.pal(n = 5, name = "PuOr"))),limits=c(-4,4)) + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  ylab('') +
  theme(axis.ticks = element_blank()) 







#Short Pretreatment vs Long pretreatment
df.R<-data.frame(unique(c(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark" & df_GSEA_specific_DESeq2$Comparison == "ShortPre_LongPre",]$pathway)))
colnames(df.R) <-"pathway"
for (i in unique(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark"& df_GSEA_specific_DESeq2$Comparison == "ShortPre_LongPre",]$Celltype)) {
  print(i)
  df.R<- merge(df.R,df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark" & df_GSEA_specific_DESeq2$Comparison == "ShortPre_LongPre" & df_GSEA_specific_DESeq2$Celltype == i,c("pathway","NES")],by="pathway", all.x = T)
}
colnames(df.R)<-c("pathway",c(unique(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark"& df_GSEA_specific_DESeq2$Comparison == "ShortPre_LongPre",]$Celltype)))
df.R[is.na(df.R)]<-0
df.R$pathway_simple <- apply(df.R,1,function(x){unlist(strsplit(x[1],"MARK_"))[2]})


df.R2<-data.frame(unique(c(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark"& df_GSEA_specific_DESeq2$Comparison == "ShortPre_LongPre",]$pathway)))
colnames(df.R2) <-"pathway"
for (i in unique(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark"& df_GSEA_specific_DESeq2$Comparison == "ShortPre_LongPre",]$Celltype)) {
  print(i)
  df.R2<- merge(df.R2,df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark" & df_GSEA_specific_DESeq2$Comparison == "ShortPre_LongPre" & df_GSEA_specific_DESeq2$Celltype == i,c("pathway","padj")],by="pathway", all.x = T)
}
colnames(df.R2)<-c("pathway",c(unique(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark"& df_GSEA_specific_DESeq2$Comparison == "ShortPre_LongPre",]$Celltype)))
df.R2[is.na(df.R2)]<-1
df.R2$pathway_simple <- apply(df.R2,1,function(x){unlist(strsplit(x[1],"MARK_"))[2]})
df.R2[,c(2:13)]<- -log(df.R2[,c(2:13)])

df.R_M <- reshape2::melt(df.R[,c(14,2:13)])
df.R2_M <- reshape2::melt(df.R2[,c(14,2:13)])
df.R3_M <- merge(df.R_M,df.R2_M, by= c("pathway_simple","variable"))

df.R3_M$variable <- factor(df.R3_M$variable, levels = c("EPC","MEP","HSC","LMPP","Memory B","CD56dim NK", "Vd1 gdT","CD4 Naive","CD4 TCM","CD4 Treg","CD8 TCM","CD8 TEMRA"))
ggplot(df.R3_M, aes(x=variable, y = pathway_simple, color = value.x, size = value.y)) + 
  geom_point() + 
  scale_color_gradientn(colours =rev(c(brewer.pal(n = 5, name = "PuOr"))),limits=c(-4,4)) + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  ylab('') +
  theme(axis.ticks = element_blank()) 




#Short Posttreatment vs Long posttreatment
df.R<-data.frame(unique(c(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark" & df_GSEA_specific_DESeq2$Comparison == "ShortPost_LongPost",]$pathway)))
colnames(df.R) <-"pathway"
for (i in unique(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark"& df_GSEA_specific_DESeq2$Comparison == "ShortPost_LongPost",]$Celltype)) {
  print(i)
  df.R<- merge(df.R,df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark" & df_GSEA_specific_DESeq2$Comparison == "ShortPost_LongPost" & df_GSEA_specific_DESeq2$Celltype == i,c("pathway","NES")],by="pathway", all.x = T)
}
colnames(df.R)<-c("pathway",c(unique(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark"& df_GSEA_specific_DESeq2$Comparison == "ShortPost_LongPost",]$Celltype)))
df.R[is.na(df.R)]<-0
df.R$pathway_simple <- apply(df.R,1,function(x){unlist(strsplit(x[1],"MARK_"))[2]})


df.R2<-data.frame(unique(c(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark"& df_GSEA_specific_DESeq2$Comparison == "ShortPost_LongPost",]$pathway)))
colnames(df.R2) <-"pathway"
for (i in unique(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark"& df_GSEA_specific_DESeq2$Comparison == "ShortPost_LongPost",]$Celltype)) {
  print(i)
  df.R2<- merge(df.R2,df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark" & df_GSEA_specific_DESeq2$Comparison == "ShortPost_LongPost" & df_GSEA_specific_DESeq2$Celltype == i,c("pathway","padj")],by="pathway", all.x = T)
}
colnames(df.R2)<-c("pathway",c(unique(df_GSEA_specific_DESeq2[df_GSEA_specific_DESeq2$gsea_type == "Hallmark"& df_GSEA_specific_DESeq2$Comparison == "ShortPost_LongPost",]$Celltype)))
df.R2[is.na(df.R2)]<-1
df.R2$pathway_simple <- apply(df.R2,1,function(x){unlist(strsplit(x[1],"MARK_"))[2]})
df.R2[,c(2:14)]<- -log(df.R2[,c(2:14)])

df.R_M <- reshape2::melt(df.R[,c(15,2:14)])
df.R2_M <- reshape2::melt(df.R2[,c(15,2:14)])
df.R3_M <- merge(df.R_M,df.R2_M, by= c("pathway_simple","variable"))

df.R3_M$variable <- factor(df.R3_M$variable, levels = c("EPC","MEP","HSC","LMPP","Memory B","CD56dim NK", "Vd1 gdT","CD4 Naive","CD4 TCM","CD4 Treg","CD8 TCM", "CD8 TEM","CD8 TEMRA"))
ggplot(df.R3_M, aes(x=variable, y = pathway_simple, color = value.x, size = value.y)) + 
  geom_point() + 
  scale_color_gradientn(colours =rev(c(brewer.pal(n = 5, name = "PuOr"))),limits=c(-4,4)) + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  ylab('') +
  theme(axis.ticks = element_blank()) 





################################### CYTOKINE EXPRESSION #########################
SASP <- c("CDKN1A","CDKN1B","CDKN2A","CDKN2B", "CDKN2C","CDKN2D","CXCL2","CXCL3","CXCL8")

DefaultAssay(BM.combined.sct.BM) = "RNA"

BM.combined.sct.BM$manualAnno_specific2 <- factor(BM.combined.sct.BM$manualAnno_specific2, levels = c("EPC","MEP","HSC","LMPP","GMP","Classical Mono","Intermediate Mono","Non-classical Mono","mDC","pDC","Naive B","Memory B","Plasmablast","Vd1 gdT","Vd2 gdT","CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL","CD4 Treg","CD8 Naive","CD8 TCM","CD8 TEM","CD8 TEMRA","MAIT","CD56dim NK","CD56bright NK",'Stromal', "Megakaryocyte","Neutrophil","Differentiating Stem Cell"))
levels(BM.combined.sct.BM$manualAnno_specific2) <- c("EPC","MEP","HSC","LMPP","GMP","Classical Mono","Intermediate Mono","Non-classical Mono","pDC","mDC","Naive B","Memory B","Plasmablast","Vd1 gdT","Vd2 gdT","CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL","CD4 Treg","CD8 Naive","CD8 TCM","CD8 TEM","CD8 TEMRA","MAIT","CD56dim NK","CD56bright NK",'Stromal', "Megakaryocyte","Neutrophil","Differentiating Stem Cell")
BM.combined.sct.BM$TimeSurvival <- factor(BM.combined.sct.BM$TimeSurvival,levels = c("Short_Pre","Short_Post","Long_Pre","Long_Post"))
Idents(BM.combined.sct.BM) <- "manualAnno_specific2"

##S100A8, S100A9
Idents(BM.combined.sct.BM) <- "manualAnno_specific2"

DefaultAssay(BM.combined.sct.BM) = "RNA"
aa<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("S100A8"),idents = c("HSC","LMPP","GMP","Classical Mono"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
bb<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("S100A9"),idents = c("HSC","LMPP","GMP","Classical Mono"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
multiplot(aa,bb,cols=1)

#CDKN & CXCL genes
aa<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("CDKN2C"),idents = c("MEP","HSC","LMPP","GMP"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
bb<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("CDKN2D"),idents = c("MEP","HSC","LMPP","GMP"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
cc<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("CXCL2"),idents = c("MEP","HSC","LMPP","GMP"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
dd<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("CXCL3"),idents = c("MEP","HSC","LMPP","GMP"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
ee<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("CXCL8"),idents = c("MEP","HSC","LMPP","GMP"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')

multiplot(cc,aa,dd,bb,ee,cols=3)

aa<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("CDKN1A"),idents = c("MEP","HSC","LMPP","GMP"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
bb<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("CDKN1B"),idents = c("MEP","HSC","LMPP","GMP"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
cc<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("CDKN2A"),idents = c("MEP","HSC","LMPP","GMP"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
dd<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("CDKN2B"),idents = c("MEP","HSC","LMPP","GMP"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
multiplot(aa,dd,bb,cc,cols=3)


#CD1 genes
DefaultAssay(BM.combined.sct.BM) = "RNA"
aa<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("CD1E"),idents = c("HSC","LMPP","GMP","Classical Mono","mDC"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
bb<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("CD1D"),idents = c("HSC","LMPP","GMP","Classical Mono","mDC"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
multiplot(aa,bb,cols=1)

DefaultAssay(BM.combined.sct.BM) = "ADT2"
aa<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("CD1d"),idents = c("mDC"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
bb<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("CD1d"),idents = c("HSC","LMPP","GMP"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')

multiplot(aa,bb,cols=1)




### SENESCENCE CDKN EXPRESSION

DefaultAssay(BM.combined.sct.BM) = "RNA"
f<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("CDKN1A"),idents = c("MEP","HSC","LMPP","GMP"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=1,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
a<-VlnPlot(object = BM.combined.sct.BM, layer = "data",  features = c("CDKN1B"),idents = c("MEP","HSC","LMPP","GMP"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=1,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
b<-VlnPlot(object = BM.combined.sct.BM, layer = "data",  features = c("CDKN2A"),idents = c("MEP","HSC","LMPP","GMP"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=1,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
c<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("CDKN2B"),idents = c("MEP","HSC","LMPP","GMP"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=1,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
d<-VlnPlot(object = BM.combined.sct.BM, layer = "data",  features = c("CDKN2C"),idents = c("MEP","HSC","LMPP","GMP"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=1,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
e<-VlnPlot(object = BM.combined.sct.BM, layer = "data",  features = c("CDKN2D"),idents = c("MEP","HSC","LMPP","GMP"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=1,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')

multiplot(f,a,b,c,d,e,cols=3)



DefaultAssay(BM.combined.sct.BM) = "RNA"
f<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("CXCL8"),idents = c("MEP","HSC","LMPP","GMP"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=1,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
a<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("CXCL2"),idents = c("MEP","HSC","LMPP","GMP"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=1,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
b<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("CXCL3"),idents = c("MEP","HSC","LMPP","GMP"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=1,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')

multiplot(a,b,f,cols=1)


#### MDSC gene signature (https://www.por-journal.com/articles/10.3389/pore.2023.1611210/full#:~:text=MDSC%2Drelated%20genes.-,We%20identified%20a%20set%20of%20six%20MDSC%2Drelated%20genes%20that,RTN4%2C%20SLC2A3%2C%20and%20TNFAIP6.)

DefaultAssay(BM.combined.sct.BM) = "RNA"
d<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("ALDOA"),idents = c("GMP","Classical Mono"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
a<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("CD52"),idents = c("GMP","Classical Mono"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
b<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("FTH1"),idents = c("GMP","Classical Mono"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
c<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("RTN4"),idents = c("GMP","Classical Mono"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
e<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("TNFAIP6"),idents = c("GMP","Classical Mono"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
f<-VlnPlot(object = BM.combined.sct.BM, layer = "data", features = c("ARG1"),idents = c("GMP","Classical Mono"), cols = c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')

multiplot(f,d,a,b,c,e,cols=2)


############### MAKE normalized Counts for PSEUDOBULKED SAMPLES TO COMPARE GENE EXPRESSION ################

DESeq2_NormCounts_broad <-data.frame()
for (broadcelltype in c("CD34")){
  print(broadcelltype)
  cluster_counts <- counts_ls_broad[[broadcelltype]]
  cluster_metadata <- metadata_ls_broad[[broadcelltype]]
  cluster_metadata$Patient <- factor(cluster_metadata$Patient)
  
  cluster_counts <- cluster_counts[,cluster_metadata$cluster_sample_id]
  
  print(all(colnames(cluster_counts) == rownames(cluster_metadata)))
  
  # Create DESeq2 object    
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ TimeSurvival)
  dds <- DESeq(dds)
  
  NormCountsData <- data.frame(counts(dds, normalize = T), stringsAsFactors = F)
  NormCountsData$Gene <- rownames(NormCountsData)
  if (broadcelltype == "CD34"){
    d <- NormCountsData
    DESeq2_NormCounts_broad <- NormCountsData }
  else {
    d2 <- NormCountsData
    DESeq2_NormCounts_broad <- merge(DESeq2_NormCounts_broad,NormCountsData,by = "Gene", all = T) }
}






DESeq2_NormCounts_specific <-data.frame()
for (specificcelltype in c("CD4 Naive","CD4 TCM","EPC","Memory B","CD8 TEMRA", "LMPP","HSC","CD56dim NK","GMP","Classical Mono","Plasmablast","mDC","pDC","CD56bright NK", "CD8 Naive","MAIT","MEP","Naive B","Vd1 gdT","Vd2 gdT","CD8 TEM","CD8 TCM","CD4 Treg")){
  tryCatch({ print(specificcelltype)
    cluster_counts <- counts_ls_specific[[specificcelltype]]
    cluster_metadata <- metadata_ls_specific[[specificcelltype]]
    cluster_metadata$Patient <- factor(cluster_metadata$Patient)
    
    cluster_counts <- cluster_counts[,cluster_metadata$cluster_sample_id]
    
    print(all(colnames(cluster_counts) == rownames(cluster_metadata)))
    
    
    # Create DESeq2 object    
    dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                  colData = cluster_metadata, 
                                  design = ~ TimeSurvival)
    dds <- DESeq(dds)
    
    NormCountsData <- data.frame(counts(dds, normalize = T), stringsAsFactors = F)
    NormCountsData$Gene <- rownames(NormCountsData)
    if (specificcelltype == "CD4 Naive"){
      d <- NormCountsData
      DESeq2_NormCounts_specific <- NormCountsData }
    else {
      d2 <- NormCountsData
      DESeq2_NormCounts_specific <- merge(DESeq2_NormCounts_specific,NormCountsData,by = "Gene", all = T) }
  })}


#################PSEUDOBULK FIGURES FOR GENE EXPRESSION###########################

######## Genes of Interest below ##########

BM.combined.sct.BM$PatientTime <- factor(BM.combined.sct.BM$PatientTime, levels = c("24780_Pre","24780_Post","122759_Pre","122759_Post","113446_Pre","113446_Post","122739_Pre","122739_Post","121535_Pre","121535_Post","121474_Pre","121474_Post","119837_Pre","119837_Post","121015_Pre","121015_Post","122976_Pre","122976_Post","114152_Pre","114152_Post","119325_Pre","119325_Post","118877_Pre","118877_Post","115585_Pre","115585_Post"))

Idents(BM.combined.sct.BM) <- "manualAnno_specific2"
GOIgenes <- c("S100A8", "S100A9")
VlnPlot(obj = BM.combined.sct.BM, features = GOIgenes, idents = c("HSC","LMPP","GMP","Classical Mono"),cols =c(mypalette,mypalette2),split.by = 'PatientTime', ncol = 1)& theme(legend.position = "top")

GOIgenes <- c("CD1D" ,"CD1E")
VlnPlot(obj = BM.combined.sct.BM, features = GOIgenes, idents = c("HSC","LMPP","GMP","Classical Mono","mDC"),cols =c(mypalette,mypalette2),split.by = 'PatientTime', ncol = 1)& theme(legend.position = "top")


GOIgenes <- c("CDKN1A","CDKN1B","CDKN2A","CDKN2B","CDKN2C","CDKN2D","CXCL2","CXCL3","CXCL8")

VlnPlot(obj = BM.combined.sct.BM, features = GOIgenes, idents = c("MEP","HSC","LMPP","GMP"),cols =c(mypalette,mypalette2),split.by = 'PatientTime', ncol = 3)& theme(legend.position = "None")




GOIgenes <- c("CDKN1A","CDKN1B","CDKN2A","CDKN2B","CDKN2C","CDKN2D","CXCL2","CXCL3","CXCL8","CLK3","ITGA5")

DESeq2_NormCounts_specific_GOI <- DESeq2_NormCounts_specific[DESeq2_NormCounts_specific$Gene %in% GOIgenes,]

m_specific_GOI <- reshape2::melt(DESeq2_NormCounts_specific_GOI, id = "Gene")
m_specific_GOI$Anno <- apply(m_specific_GOI, 1, function(x) {unlist(strsplit(x[2], "_"))[1]})
m_specific_GOI$Patient <- apply(m_specific_GOI, 1, function(x) {unlist(strsplit(x[2], "_"))[2]})
m_specific_GOI$Survival <- apply(m_specific_GOI, 1, function(x) {unlist(strsplit(x[2], "_"))[3]})
m_specific_GOI$Time <- apply(m_specific_GOI, 1, function(x) {unlist(strsplit(x[2], "_"))[4]})
m_specific_GOI$SurvivalTime <- paste(m_specific_GOI$Survival,m_specific_GOI$Time, sep = "-")
m_specific_GOI$PatientSurvivalTime <- paste(m_specific_GOI$Patient,m_specific_GOI$Survival,m_specific_GOI$Time, sep = "-")


m_specific_GOI$Anno <- factor(m_specific_GOI$Anno, levels = c("EPC","MEP","HSC","LMPP","GMP","Classical.Mono","Intermediate.Mono","Non-classical.Mono","mDC","pDC","Naive.B","Memory.B","Plasmablast","Vd1.gdT","Vd2.gdT","CD4.Naive","CD4.TCM","CD4.TEM","CD4.CTL","CD4.Treg","CD8.Naive","CD8.TCM","CD8.TEM","CD8.TEMRA","MAIT","CD56dim.NK","CD56bright.NK",'Stromal',"Megakaryocyte","Neutrophil","Differentiating.Stem.Cell"))
m_specific_GOI$SurvivalTime <- factor(m_specific_GOI$SurvivalTime, levels = c("Short-Pre","Short-Post","Long-Pre","Long-Post"))

m_specific_GOI$Gene <- factor(m_specific_GOI$Gene, levels = GOIgenes)


p<-ggplot(m_specific_GOI[m_specific_GOI$Anno %in% c("MEP","HSC","LMPP","GMP"),])+geom_boxplot(aes(Anno,value, color=SurvivalTime),outlier.shape = NA,width = 0.8)+geom_point(aes(Anno,value, group = SurvivalTime, color=SurvivalTime),size=1.5,position = position_dodge(width=0.8))+ggtitle("Immune Frequency Dynamics (Broad)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Anno", y = "normCounts")+scale_fill_manual(values =c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"))+scale_color_manual(values =c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
p+ facet_wrap( ~ Gene, scales="free",ncol=6)

pp<-ggplot(m_specific_GOI[m_specific_GOI$Anno %in% c("HSC","LMPP","GMP","Classical.Mono","mDC"),])+geom_boxplot(aes(Anno,value, color=SurvivalTime),outlier.shape = NA,width = 0.4)+geom_point(aes(Anno,value, group = SurvivalTime, color=Patient),size=3,position = position_dodge(width=0.4))+ggtitle("Immune Frequency Dynamics (Broad)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Anno", y = "normCounts")+scale_fill_manual(values =c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3",mypalette,mypalette2))+scale_color_manual(values =c("#d4d0d1","#5f5e60","#e6c6e4","#8247a3",mypalette,mypalette2))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
pp+ facet_wrap( ~ Gene, scales="free",ncol=2)





############## LOOK AT TCR/BCR REPRETOIRE USING TRUST4 ############

BM.combined.sct <- LoadSeuratRds("BM.combined.sct_2025.rds")

library(DescTools)

BM.combined.sct$SampleInfo_Type_Time <- paste(BM.combined.sct$Patient,BM.combined.sct$SampleType,BM.combined.sct$Time,sep="-")
BM.combined.sct$SampleInfo_Type_Time<- factor(BM.combined.sct$SampleInfo_Type_Time, levels = c("24780-BM-Scr","24780-BM-PC2","122759-BM-Scr","122759-BM-PC2","113446-BM-Scr","113446-BM-PC2","122739-BM-Scr","122739-BM-PC2","121535-BM-Scr","121535-BM-PC2","121474-BM-Scr","121474-BM-PC2","119837-BM-Scr","119837-BM-PC2","121015-BM-Scr","121015-BM-PC2","122976-BM-Scr","122976-BM-PC2","114152-BM-Scr","114152-BM-PC2","119325-BM-Scr","119325-BM-PC2","118877-BM-Scr","118877-BM-PC2","115585-BM-Scr","115585-BM-PC2","24780-PBMC-Scr","24780-PBMC-C2D1","24780-PBMC-C2D8","121474-PBMC-Scr","121474-PBMC-C1D8","121474-PBMC-C2D1","121015-PBMC-Scr","121015-PBMC-C1D8","121015-PBMC-C2D1","118877-PBMC-Scr","118877-PBMC-C1D8","118877-PBMC-C2D1"))

BM.combined.sct@meta.data$PatientTime <- paste(BM.combined.sct@meta.data$Patient,BM.combined.sct@meta.data$Time2, sep="_")
BM.combined.sct@meta.data$TimeSurvival <- paste(BM.combined.sct@meta.data$Survival,BM.combined.sct@meta.data$Time2,sep = "_")

BM.combined.sct$SampleInfo_Type_Time_Specific <- paste(BM.combined.sct$manualAnno_specific2,BM.combined.sct$SampleInfo_Type_Time, sep ="_")
BM.combined.sct$SampleInfo_Type_Time_broad <- paste(BM.combined.sct$manualAnno_specific,BM.combined.sct$SampleInfo_Type_Time, sep ="_")

#Filtered TCR/BCR in "MDS_Exp_barcode_report.tsv" to identify major chains in each cell.

setwd("./MDS_TRUST4/vdjtools")

exp2_report <- read.table("MDS_Exp2_barcode_report.tsv", sep = "\t", header = F, stringsAsFactors = F)
exp3_report <- read.table("MDS_Exp3_barcode_report.tsv", sep = "\t", header = F, stringsAsFactors = F)
exp4_report <- read.table("MDS_Exp4_barcode_report.tsv", sep = "\t", header = F, stringsAsFactors = F)
exp5_report <- read.table("MDS_Exp5_barcode_report.tsv", sep = "\t", header = F, stringsAsFactors = F)
exp6_report <- read.table("MDS_Exp6_barcode_report.tsv", sep = "\t", header = F, stringsAsFactors = F)
exp7_report <- read.table("MDS_Exp7_barcode_report.tsv", sep = "\t", header = F, stringsAsFactors = F)
exp7b_report <- read.table("MDS_Exp7b_barcode_report.tsv", sep = "\t", header = F, stringsAsFactors = F)
exp8_report <- read.table("MDS_Exp8_barcode_report.tsv", sep = "\t", header = F, stringsAsFactors = F)

exp2_report$CellID <- paste("Exp2", exp2_report$V1, sep = "_")
exp3_report$CellID <- paste("Exp3", exp3_report$V1, sep = "_")
exp4_report$CellID <- paste("Exp4", exp4_report$V1, sep = "_")
exp5_report$CellID <- paste("Exp5", exp5_report$V1, sep = "_")
exp6_report$CellID <- paste("Exp6", exp6_report$V1, sep = "_")
exp7_report$CellID <- paste("Exp7", exp7_report$V1, sep = "_")
exp7b_report$CellID <- paste("Exp7b", exp7b_report$V1, sep = "_")
exp8_report$CellID <- paste("Exp8", exp8_report$V1, sep = "_")

MDS_allExp_report <- rbind(exp2_report[,c(7,2:6)],exp3_report[,c(7,2:6)],exp4_report[,c(7,2:6)],exp5_report[,c(7,2:6)],exp6_report[,c(7,2:6)],exp7_report[,c(7,2:6)],exp7b_report[,c(7,2:6)],exp8_report[,c(7,2:6)])
colnames(MDS_allExp_report) <- c("barcode",	"cell_type",	"chain1",	"chain2",	"secondary_chain1",	"secondary_chain2")




#more data on cell TCR/BCR calls
setwd("./MDS_TRUST4/airr")
exp2 <- read.table("MDS_Exp2_barcode_airr.tsv", sep = "\t", header = T, stringsAsFactors = F)
exp3 <- read.table("MDS_Exp3_barcode_airr.tsv", sep = "\t", header = T, stringsAsFactors = F)
exp4 <- read.table("MDS_Exp4_barcode_airr.tsv", sep = "\t", header = T, stringsAsFactors = F)
exp5 <- read.table("MDS_Exp5_barcode_airr.tsv", sep = "\t", header = T, stringsAsFactors = F)
exp6 <- read.table("MDS_Exp6_barcode_airr.tsv", sep = "\t", header = T, stringsAsFactors = F)
exp7 <- read.table("MDS_Exp7_barcode_airr.tsv", sep = "\t", header = T, stringsAsFactors = F)
exp7b <- read.table("MDS_Exp7b_barcode_airr.tsv", sep = "\t", header = T, stringsAsFactors = F)
exp8 <- read.table("MDS_Exp8_barcode_airr.tsv", sep = "\t", header = T, stringsAsFactors = F)

exp2$CellID <- paste("Exp2", exp2$cell_id, sep = "_")
exp3$CellID <- paste("Exp3", exp3$cell_id, sep = "_")
exp4$CellID <- paste("Exp4", exp4$cell_id, sep = "_")
exp5$CellID <- paste("Exp5", exp5$cell_id, sep = "_")
exp6$CellID <- paste("Exp6", exp6$cell_id, sep = "_")
exp7$CellID <- paste("Exp7", exp7$cell_id, sep = "_")
exp7b$CellID <- paste("Exp7b", exp7b$cell_id, sep = "_")
exp8$CellID <- paste("Exp8", exp8$cell_id, sep = "_")



MDS_allExp_airr <- rbind(exp2[,c(24,2:23)],exp3[,c(24,2:23)],exp4[,c(24,2:23)],exp5[,c(24,2:23)],exp6[,c(24,2:23)],exp7[,c(24,2:23)],exp7b[,c(24,2:23)],exp8[,c(24,2:23)])

MDS_allExp_airr[MDS_allExp_airr$complete_vdj == "F",]

#Filter cells with complete VDJ reconstruction

MDS_allExp_report_complete <- MDS_allExp_report

MDS_allExp_report_complete_abT <- MDS_allExp_report_complete[MDS_allExp_report_complete$cell_type == "abT",]
MDS_allExp_report_complete_B <- MDS_allExp_report_complete[MDS_allExp_report_complete$cell_type == "B",]
MDS_allExp_report_complete_gdT <- MDS_allExp_report_complete[MDS_allExp_report_complete$cell_type == "gdT",]


MDS_allExp_report_complete_abT$chain1a <- apply(MDS_allExp_report_complete_abT, 1, function(x) {unlist(strsplit(x[3],","))[1]})
MDS_allExp_report_complete_abT$chain1b <- apply(MDS_allExp_report_complete_abT, 1, function(x) {unlist(strsplit(x[3],","))[2]})
MDS_allExp_report_complete_abT$chain1c <- apply(MDS_allExp_report_complete_abT, 1, function(x) {unlist(strsplit(x[3],","))[3]})
MDS_allExp_report_complete_abT$chain1d <- apply(MDS_allExp_report_complete_abT, 1, function(x) {unlist(strsplit(x[3],","))[4]})
MDS_allExp_report_complete_abT$chain1e <- apply(MDS_allExp_report_complete_abT, 1, function(x) {unlist(strsplit(x[3],","))[5]})
MDS_allExp_report_complete_abT$chain1f <- apply(MDS_allExp_report_complete_abT, 1, function(x) {unlist(strsplit(x[3],","))[6]})


MDS_allExp_report_complete_abT$chain2a <- apply(MDS_allExp_report_complete_abT, 1, function(x) {unlist(strsplit(x[4],","))[1]})
MDS_allExp_report_complete_abT$chain2b <- apply(MDS_allExp_report_complete_abT, 1, function(x) {unlist(strsplit(x[4],","))[2]})
MDS_allExp_report_complete_abT$chain2c <- apply(MDS_allExp_report_complete_abT, 1, function(x) {unlist(strsplit(x[4],","))[3]})
MDS_allExp_report_complete_abT$chain2d <- apply(MDS_allExp_report_complete_abT, 1, function(x) {unlist(strsplit(x[4],","))[4]})
MDS_allExp_report_complete_abT$chain2e <- apply(MDS_allExp_report_complete_abT, 1, function(x) {unlist(strsplit(x[4],","))[5]})
MDS_allExp_report_complete_abT$chain2f <- apply(MDS_allExp_report_complete_abT, 1, function(x) {unlist(strsplit(x[4],","))[6]})


MDS_allExp_report_complete_abT$CDR3aa <- paste(MDS_allExp_report_complete_abT$chain1f,MDS_allExp_report_complete_abT$chain2f,sep = "+")

t3 <-table(MDS_allExp_report_complete_abT$CDR3aa)
t3[order(-t3)]






MDS_allExp_report_complete_B$chain1a <- apply(MDS_allExp_report_complete_B, 1, function(x) {unlist(strsplit(x[3],","))[1]})
MDS_allExp_report_complete_B$chain1b <- apply(MDS_allExp_report_complete_B, 1, function(x) {unlist(strsplit(x[3],","))[2]})
MDS_allExp_report_complete_B$chain1c <- apply(MDS_allExp_report_complete_B, 1, function(x) {unlist(strsplit(x[3],","))[3]})
MDS_allExp_report_complete_B$chain1d <- apply(MDS_allExp_report_complete_B, 1, function(x) {unlist(strsplit(x[3],","))[4]})
MDS_allExp_report_complete_B$chain1e <- apply(MDS_allExp_report_complete_B, 1, function(x) {unlist(strsplit(x[3],","))[5]})
MDS_allExp_report_complete_B$chain1f <- apply(MDS_allExp_report_complete_B, 1, function(x) {unlist(strsplit(x[3],","))[6]})

MDS_allExp_report_complete_B$chain2a <- apply(MDS_allExp_report_complete_B, 1, function(x) {unlist(strsplit(x[4],","))[1]})
MDS_allExp_report_complete_B$chain2b <- apply(MDS_allExp_report_complete_B, 1, function(x) {unlist(strsplit(x[4],","))[2]})
MDS_allExp_report_complete_B$chain2c <- apply(MDS_allExp_report_complete_B, 1, function(x) {unlist(strsplit(x[4],","))[3]})
MDS_allExp_report_complete_B$chain2d <- apply(MDS_allExp_report_complete_B, 1, function(x) {unlist(strsplit(x[4],","))[4]})
MDS_allExp_report_complete_B$chain2e <- apply(MDS_allExp_report_complete_B, 1, function(x) {unlist(strsplit(x[4],","))[5]})
MDS_allExp_report_complete_B$chain2f <- apply(MDS_allExp_report_complete_B, 1, function(x) {unlist(strsplit(x[4],","))[6]})

MDS_allExp_report_complete_B$CDR3aa <- paste(MDS_allExp_report_complete_B$chain1f,MDS_allExp_report_complete_B$chain2f,sep = "+")

t3 <-table(MDS_allExp_report_complete_B$CDR3aa)
t3[order(-t3)]




MDS_allExp_report_complete_gdT$chain1a <- apply(MDS_allExp_report_complete_gdT, 1, function(x) {unlist(strsplit(x[3],","))[1]})
MDS_allExp_report_complete_gdT$chain1b <- apply(MDS_allExp_report_complete_gdT, 1, function(x) {unlist(strsplit(x[3],","))[2]})
MDS_allExp_report_complete_gdT$chain1c <- apply(MDS_allExp_report_complete_gdT, 1, function(x) {unlist(strsplit(x[3],","))[3]})
MDS_allExp_report_complete_gdT$chain1d <- apply(MDS_allExp_report_complete_gdT, 1, function(x) {unlist(strsplit(x[3],","))[4]})
MDS_allExp_report_complete_gdT$chain1e <- apply(MDS_allExp_report_complete_gdT, 1, function(x) {unlist(strsplit(x[3],","))[5]})
MDS_allExp_report_complete_gdT$chain1f <- apply(MDS_allExp_report_complete_gdT, 1, function(x) {unlist(strsplit(x[3],","))[6]})

MDS_allExp_report_complete_gdT$chain2a <- apply(MDS_allExp_report_complete_gdT, 1, function(x) {unlist(strsplit(x[4],","))[1]})
MDS_allExp_report_complete_gdT$chain2b <- apply(MDS_allExp_report_complete_gdT, 1, function(x) {unlist(strsplit(x[4],","))[2]})
MDS_allExp_report_complete_gdT$chain2c <- apply(MDS_allExp_report_complete_gdT, 1, function(x) {unlist(strsplit(x[4],","))[3]})
MDS_allExp_report_complete_gdT$chain2d <- apply(MDS_allExp_report_complete_gdT, 1, function(x) {unlist(strsplit(x[4],","))[4]})
MDS_allExp_report_complete_gdT$chain2e <- apply(MDS_allExp_report_complete_gdT, 1, function(x) {unlist(strsplit(x[4],","))[5]})
MDS_allExp_report_complete_gdT$chain2f <- apply(MDS_allExp_report_complete_gdT, 1, function(x) {unlist(strsplit(x[4],","))[6]})

MDS_allExp_report_complete_gdT$CDR3aa <- paste(MDS_allExp_report_complete_gdT$chain1f,MDS_allExp_report_complete_gdT$chain2f,sep = "+")


t3 <-table(MDS_allExp_report_complete_gdT$CDR3aa)
t3[order(-t3)]


MDS_allExp_report_complete2 <- rbind(MDS_allExp_report_complete_abT[,c("barcode","cell_type","CDR3aa")],MDS_allExp_report_complete_B[,c("barcode","cell_type","CDR3aa")],MDS_allExp_report_complete_gdT[,c("barcode","cell_type","CDR3aa")])
colnames(MDS_allExp_report_complete2)[1] <- "V1"

MDS_allExp_report_complete2$CDR3aa <- apply(MDS_allExp_report_complete2,1, function(x) {
  if (grepl("\\+NA", x[3]) == "TRUE") {unlist(strsplit(x[3],"\\+"))[1]} 
  else if (grepl("^NA\\+", x[3]) == "TRUE") {
    unlist(strsplit(x[3],"\\+"))[2]} 
  else {x[3]}
})

MDS_allExp_report_complete2$ChainDetect <- apply(MDS_allExp_report_complete2,1, function(x) {
  if (grepl("\\+", x[3]) == "TRUE") {"Double"} 
  else {"Single"} 
})


#MERGE THIS INFO WITH CLINICAL METADATA
library(dplyr)
library(tidyr)

met <- read.table("BMmerged_metadata_HTOclassification_mito15_wClinical_Celltype.txt", sep = "\t", header = F, stringsAsFactors = F)

met1 <- merge(met,MDS_allExp_report_complete2, by = "V1", all.x = T)
met1$cell_type <- met1$cell_type %>% replace_na('NotDetected')
met1$CDR3aa <- met1$CDR3aa %>% replace_na('NotDetected')
met1$ChainDetect <- met1$ChainDetect %>% replace_na('NotDetected')


#filter out discordant read calls for cell type between Seurat and TRUST4
met1[met1$cell_type == "abT" & !(met1$V47 %in% c("CD4","CD8")),]$CDR3aa <- "NotDetected"
met1[met1$cell_type == "abT" & !(met1$V47 %in% c("CD4","CD8")),]$ChainDetect <- "NotDetected"

met1[met1$cell_type == "gdT" & !(met1$V47 %in% c("gdT")),]$CDR3aa <- "NotDetected"
met1[met1$cell_type == "gdT" & !(met1$V47 %in% c("gdt")),]$ChainDetect <- "NotDetected"

met1[met1$cell_type == "B" & !(met1$V47 %in% c("B cell", "Plasmablast")),]$CDR3aa <- "NotDetected"
met1[met1$cell_type == "B" & !(met1$V47 %in% c("B cell", "Plasmablast")),]$ChainDetect <- "NotDetected"

met1$broad3 <- met1$V47
met1[met1$broad3 == "Plasmablast",]$broad3 <- "B cell" #combine B cell and Plasmablast for BCR analysis

setwd("/varidata/research/projects/jang/JJANG/projects/CITEseq/VAI_OnDemand_Rstudio/Seurat5_nolayers")
write.table(met1, "BMmerged_metadata_HTOclassification_mito15_wClinical_121823_UpdatedTRUST4_012224.txt", col.names = F, row.names = F, quote = F, sep = "\t")




####
metaclinical <- read.table("BMmerged_metadata_HTOclassification_mito15_wClinical_UpdatedTRUST4.txt", sep = "\t", header = F, stringsAsFactors = F)
rownames(metaclinical) <- metaclinical$V1
metaclinical <- metaclinical[,c(2:58)]
colnames(metaclinical) <- c(colnames(BM.combined.sct@meta.data),"TRUST4_CellType","ChainInfo","ChainCount","Detection","broad3")

BM.combined.sct <-AddMetaData(object = BM.combined.sct, metadata = metaclinical)
BM.combined.sct$SampleInfo_Type_Time_broad <- paste(BM.combined.sct$broad3,BM.combined.sct$SampleInfo_Type_Time, sep ="_") #replace plasmablast with B cell for BCR analysis


BM.combined.sct.BM <- subset(BM.combined.sct, subset = SampleType == "BM")
BM.combined.sct.PBMC <- subset(BM.combined.sct, subset = SampleType == "PBMC")


######## Filter and aggregate TCR/BCR info for bone marrow cells ####

CD4_TCR <- metaclinical[metaclinical$SampleType == "BM" & metaclinical$broad3 == "CD4" & metaclinical$ChainInfo == "abT",c("Patient","Time2","Survival","PatientTime","manualAnno_broad2" ,"SampleInfo_Type_Time","TimeSurvival","ChainInfo","ChainCount","broad3")]

CD4_TCR_summary = data.frame()
for (i in unique(CD4_TCR$PatientTime)) {
  d <- data.frame(table(CD4_TCR[CD4_TCR$PatientTime == i,]$ChainCount))
  d$total <- sum(d$Freq)
  d$CellType <- "CD4"
  d$PatientTime <- i
  d$Patient <- unique(CD4_TCR[CD4_TCR$PatientTime == i,]$Patient)
  d$Time2 <- unique(CD4_TCR[CD4_TCR$PatientTime == i,]$Time2)
  d$Survival<- unique(CD4_TCR[CD4_TCR$PatientTime == i,]$Survival)
  d$TimeSurvival <- unique(CD4_TCR[CD4_TCR$PatientTime == i,]$TimeSurvival)
  d$SampleInfo_Type_Time<- unique(CD4_TCR[CD4_TCR$PatientTime == i,]$SampleInfo_Type_Time)
  d$broad3<- unique(CD4_TCR[CD4_TCR$PatientTime == i,]$broad3)
  CD4_TCR_summary <- rbind(CD4_TCR_summary,d)
}

CD8_TCR <- metaclinical[metaclinical$SampleType == "BM" & metaclinical$broad3 == "CD8" & metaclinical$ChainInfo == "abT",c("Patient","Time2","Survival","PatientTime","manualAnno_broad2","SampleInfo_Type_Time","TimeSurvival","ChainInfo","ChainCount","broad3")]

CD8_TCR_summary = data.frame()
for (i in unique(CD8_TCR$PatientTime)) {
  d <- data.frame(table(CD8_TCR[CD8_TCR$PatientTime == i,]$ChainCount))
  d$total <- sum(d$Freq)
  d$CellType <- "CD8"
  d$PatientTime <- i
  d$Patient <- unique(CD8_TCR[CD8_TCR$PatientTime == i,]$Patient)
  d$Time2 <- unique(CD8_TCR[CD8_TCR$PatientTime == i,]$Time2)
  d$Survival<- unique(CD8_TCR[CD8_TCR$PatientTime == i,]$Survival)
  d$TimeSurvival <- unique(CD8_TCR[CD8_TCR$PatientTime == i,]$TimeSurvival)
  d$SampleInfo_Type_Time<- unique(CD8_TCR[CD8_TCR$PatientTime == i,]$SampleInfo_Type_Time)
  d$broad3<- unique(CD8_TCR[CD8_TCR$PatientTime == i,]$broad3)
  CD8_TCR_summary <- rbind(CD8_TCR_summary,d)
}


B_BCR <- metaclinical[metaclinical$SampleType == "BM" & metaclinical$broad3 == "B cell" & metaclinical$ChainInfo == "B",c("Patient","Time2","Survival","PatientTime","manualAnno_broad2" ,"SampleInfo_Type_Time","TimeSurvival","ChainInfo","ChainCount","broad3")]

B_BCR_summary = data.frame()
for (i in unique(B_BCR$PatientTime)) {
  d <- data.frame(table(B_BCR[B_BCR$PatientTime == i,]$ChainCount))
  d$total <- sum(d$Freq)
  d$CellType <- "B cell"
  d$PatientTime <- i
  d$Patient <- unique(B_BCR[B_BCR$PatientTime == i,]$Patient)
  d$Time2 <- unique(B_BCR[B_BCR$PatientTime == i,]$Time2)
  d$Survival<- unique(B_BCR[B_BCR$PatientTime == i,]$Survival)
  d$TimeSurvival <- unique(B_BCR[B_BCR$PatientTime == i,]$TimeSurvival)
  d$SampleInfo_Type_Time<- unique(B_BCR[B_BCR$PatientTime == i,]$SampleInfo_Type_Time)
  d$broad3<- unique(B_BCR[B_BCR$PatientTime == i,]$broad3)
  B_BCR_summary <- rbind(B_BCR_summary,d)
}


gdT_TCR <- metaclinical[metaclinical$SampleType == "BM" & metaclinical$broad3 == "gdT" & metaclinical$ChainInfo == "gdT",c("Patient","Time2","Survival","PatientTime","manualAnno_broad2","SampleInfo_Type_Time","TimeSurvival","ChainInfo","ChainCount","broad3")]

gdT_TCR_summary = data.frame()
for (i in unique(gdT_TCR$PatientTime)) {
  d <- data.frame(table(gdT_TCR[gdT_TCR$PatientTime == i,]$ChainCount))
  d$total <- sum(d$Freq)
  d$CellType <- "gdT"
  d$PatientTime <- i
  d$Patient <- unique(gdT_TCR[gdT_TCR$PatientTime == i,]$Patient)
  d$Time2 <- unique(gdT_TCR[gdT_TCR$PatientTime == i,]$Time2)
  d$Survival<- unique(gdT_TCR[gdT_TCR$PatientTime == i,]$Survival)
  d$TimeSurvival <- unique(gdT_TCR[gdT_TCR$PatientTime == i,]$TimeSurvival)
  d$SampleInfo_Type_Time<- unique(gdT_TCR[gdT_TCR$PatientTime == i,]$SampleInfo_Type_Time)
  d$broad3<- unique(gdT_TCR[gdT_TCR$PatientTime == i,]$broad3)
  gdT_TCR_summary <- rbind(gdT_TCR_summary,d)
}



#Combine all data
t3 <- rbind(gdT_TCR_summary,B_BCR_summary,CD4_TCR_summary,CD8_TCR_summary)

t3$CellType <- factor(t3$CellType, levels = c("gdT","CD4","CD8","B cell"))

t3$Patient <- factor(t3$Patient, levels = c("24780","122759","113446","122739","121535","121474","119837","121015","122976","114152","119325","118877","115585"))

t3$SampleInfo_Type_Time <-factor(t3$SampleInfo_Type_Time,levels = c("24780-BM-Scr","24780-BM-PC2","122759-BM-Scr","122759-BM-PC2","113446-BM-Scr","113446-BM-PC2","122739-BM-Scr","122739-BM-PC2","121535-BM-Scr","121535-BM-PC2","121474-BM-Scr","121474-BM-PC2","119837-BM-Scr","119837-BM-PC2","121015-BM-Scr","121015-BM-PC2","122976-BM-Scr","122976-BM-PC2","114152-BM-Scr","114152-BM-PC2","119325-BM-Scr","119325-BM-PC2","118877-BM-Scr","118877-BM-PC2","115585-BM-Scr","115585-BM-PC2"))
t3$Time2 <- factor(t3$Time2, levels = c("Pre","Post"))
t3$Survival <- factor(t3$Survival, levels = c("Short","Long"))
t3$TimeSurvival <- factor(t3$TimeSurvival, levels = c("Short_Pre","Short_Post","Long_Pre","Long_Post"))
t3<- t3[order(t3$Patient),]

t3_total <- unique(t3[,c("total","SampleInfo_Type_Time","CellType")])

p<-ggplot(t3_total, aes(x=SampleInfo_Type_Time, y=total, fill=CellType)) +geom_bar(stat="identity",colour = "black", position="fill")+ggtitle("TCR/BCR chain Detection")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Patient Samples", y = "Frequency")+scale_fill_manual(values = c('#ff9966',"#339900",'#bae1ff','#c27ba0'))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))

pp<-ggplot(t3_total, aes(x=SampleInfo_Type_Time, y=total, fill=CellType)) +geom_bar(stat="identity",colour = "black", position="stack")+ggtitle("TCR/BCR chain Detection")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Patient Samples", y = "Frequency")+scale_fill_manual(values = c('#ff9966',"#339900",'#bae1ff','#c27ba0'))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))

multiplot(pp,p,cols = 1)


#Shannon's Entropy (load immunarch for "entropy" function)
library(immunarch)
shannon_data <- data.frame()

for (i in unique(t3$broad3)) {
  for (x in unique(t3[t3$broad3 == i,]$SampleInfo_Type_Time)) {
    sE <- entropy(t3[t3$broad3 == i & t3$SampleInfo_Type_Time == x,c("Freq")],.do.norm = NA)
    
    d <- data.frame(c(i),c(x),c(sE),c(unique(t3[t3$broad3 == i & t3$SampleInfo_Type_Time == x,]$Time2)),c(unique(t3[t3$broad3 == i & t3$SampleInfo_Type_Time == x,]$Survival)),c(unique(t3[t3$broad3 == i & t3$SampleInfo_Type_Time == x,]$Patient)),c(unique(t3[t3$broad3 == i & t3$SampleInfo_Type_Time == x,]$TimeSurvival)))
    shannon_data <- rbind(shannon_data,d)
    } 
}
colnames(shannon_data) <- c("CellType","SampleInfo_Type_Time","shannonscore","Time2","Survival","Patient","TimeSurvival")



PatientColor <- c("#111111","#555555","#999999","#f1ddf1","#eaccea","#e3bbe3","#ddaadd","#c699c6","#de94f5","#cb53ef","#be29ec","#9820bc","#851ca5")
p<-ggplot(shannon_data)+geom_boxplot(aes(Time2,shannonscore, group = interaction(Time2, Survival), color=Survival),outlier.shape = NA,width = 0.4)+geom_line(aes(Time2,shannonscore, group = Patient, color=Patient),lwd=1,alpha=0.3)+geom_point(aes(Time2,shannonscore,group=interaction(Time2, Survival), color=factor(Patient)),size=3,position = position_dodge(width=0.4))+ggtitle("Immune Chain (Shannon Entropy)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Treatment", y = "Score")+scale_fill_manual(values =c("#555555","#c699c6",PatientColor))+scale_color_manual(values =c("#555555","#c699c6",PatientColor))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
p+ facet_wrap( ~ CellType, scales="free",ncol=4)

## Not enough gdT cells for robust analysis



######## Check for clonal expansion for top clones (>3 clones ) in bone marrow #################

t3_topclones <- t3[t3$Freq > 3,]

t3_topclones_CD4 <- t3_topclones[t3_topclones$CellType == "CD4",]

t3_topclones_CD4all <- t3[t3$CellType == "CD4" & t3$Var1 %in% t3_topclones_CD4$Var1,]
t3_topclones_CD4all$percent <- t3_topclones_CD4all$Freq/t3_topclones_CD4all$total

CD4_clonal <- dcast(data.table(t3_topclones_CD4all[,c("Var1","SampleInfo_Type_Time","percent")]), Var1 ~ SampleInfo_Type_Time,value.var = "percent")
CD4_clonal1 <- CD4_clonal[,-1]
CD4_patient <- unique(unlist(strsplit(colnames(CD4_clonal1),'-'))[3*(1:length(colnames(CD4_clonal1)))-2])

CD4_final <- t3_topclones_CD4all[t3_topclones_CD4all$Patient %in% CD4_patient,c("Var1","Patient","Time2","Survival","percent")]
for (i in unique(CD4_final$Patient)) {
  for (x in CD4_final[CD4_final$Patient == i,]$Var1)
    if (length(CD4_final[CD4_final$Patient == i & CD4_final$Var1 == x,]$Time2) <2) {
      if (CD4_final[CD4_final$Patient == i & CD4_final$Var1 == x,]$Time2 == "Pre") {
        newdata <- data.frame(x,i,"Post",CD4_final[CD4_final$Patient == i & CD4_final$Var1 == x,]$Survival,0) 
      } else {newdata <- data.frame(x,i,"Pre",CD4_final[CD4_final$Patient == i & CD4_final$Var1 == x,]$Survival,0) }
      colnames(newdata) <- colnames(CD4_final)
      CD4_final <- rbind(CD4_final, newdata)  
    } 
}

p0 <- ggplot(CD4_final)+geom_line(aes(Time2,percent, group  = interaction(Var1,Patient), color=Survival),lwd=1,alpha=0.3)+geom_point(aes(Time2,percent,group=interaction(Time2, Survival), color=Survival),size=3,position = position_dodge(width=0.4))+ggtitle("CD4 Immune Chain (clonal)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Treatment", y = "Clonal fraction")+scale_fill_manual(values =c("#555555","#c699c6"))+scale_color_manual(values =c("#555555","#c699c6"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
p0 + facet_wrap( ~ Patient, scales="free",ncol=4)





##CD8
t3_topclones_CD8 <- t3_topclones[t3_topclones$CellType == "CD8",]

t3_topclones_CD8all <- t3[t3$CellType == "CD8" & t3$Var1 %in% t3_topclones_CD8$Var1,]
t3_topclones_CD8all$percent <- t3_topclones_CD8all$Freq/t3_topclones_CD8all$total

CD8_clonal <- dcast(data.table(t3_topclones_CD8all[,c("Var1","SampleInfo_Type_Time","percent")]), Var1 ~ SampleInfo_Type_Time,value.var = "percent")
CD8_clonal1 <- CD8_clonal[,-1]
CD8_patient <- unique(unlist(strsplit(colnames(CD8_clonal1),'-'))[3*(1:length(colnames(CD8_clonal1)))-2])

CD8_final <- t3_topclones_CD8all[t3_topclones_CD8all$Patient %in% CD8_patient,c("Var1","Patient","Time2","Survival","percent")]
for (i in unique(CD8_final$Patient)) {
  for (x in CD8_final[CD8_final$Patient == i,]$Var1)
    if (length(CD8_final[CD8_final$Patient == i & CD8_final$Var1 == x,]$Time2) <2) {
      if (CD8_final[CD8_final$Patient == i & CD8_final$Var1 == x,]$Time2 == "Pre") {
        newdata <- data.frame(x,i,"Post",CD8_final[CD8_final$Patient == i & CD8_final$Var1 == x,]$Survival,0) 
      } else {newdata <- data.frame(x,i,"Pre",CD8_final[CD8_final$Patient == i & CD8_final$Var1 == x,]$Survival,0) }
      colnames(newdata) <- colnames(CD8_final)
      CD8_final <- rbind(CD8_final, newdata)  
    } 
}

p1 <- ggplot(CD8_final)+geom_line(aes(Time2,percent, group  = interaction(Var1,Patient), color=Survival),lwd=1,alpha=0.3)+geom_point(aes(Time2,percent,group=interaction(Time2, Survival), color=Survival),size=3,position = position_dodge(width=0.4))+ggtitle("CD8 Immune Chain (clonal)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Treatment", y = "Clonal fraction")+scale_fill_manual(values =c("#555555","#c699c6"))+scale_color_manual(values =c("#555555","#c699c6"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
p1 + facet_wrap( ~ Patient, scales="free",ncol=4)



t3_topclones_B <- t3_topclones[t3_topclones$CellType == "B cell",]

t3_topclones_Ball <- t3[t3$CellType == "B cell" & t3$Var1 %in% t3_topclones_B$Var1,]
t3_topclones_Ball$percent <- t3_topclones_Ball$Freq/t3_topclones_Ball$total

B_clonal <- dcast(data.table(t3_topclones_Ball[,c("Var1","SampleInfo_Type_Time","percent")]), Var1 ~ SampleInfo_Type_Time,value.var = "percent")
B_clonal1 <- B_clonal[,-1]
B_patient <- unique(unlist(strsplit(colnames(B_clonal1),'-'))[3*(1:length(colnames(B_clonal1)))-2])

B_final <- t3_topclones_Ball[t3_topclones_Ball$Patient %in% B_patient,c("Var1","Patient","Time2","Survival","percent")]
for (i in unique(B_final$Patient)) {
  for (x in B_final[B_final$Patient == i,]$Var1)
    if (length(B_final[B_final$Patient == i & B_final$Var1 == x,]$Time2) <2) {
      if (B_final[B_final$Patient == i & B_final$Var1 == x,]$Time2 == "Pre") {
        newdata <- data.frame(x,i,"Post",B_final[B_final$Patient == i & B_final$Var1 == x,]$Survival,0) 
      } else {newdata <- data.frame(x,i,"Pre",B_final[B_final$Patient == i & B_final$Var1 == x,]$Survival,0) }
      colnames(newdata) <- colnames(B_final)
      B_final <- rbind(B_final, newdata)  
    } 
}

p2 <- ggplot(B_final)+geom_line(aes(Time2,percent, group  = interaction(Var1,Patient), color=Survival),lwd=1,alpha=0.3)+geom_point(aes(Time2,percent,group=interaction(Time2, Survival), color=Survival),size=3,position = position_dodge(width=0.4))+ggtitle("B Immune Chain (clonal)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Treatment", y = "Clonal fraction")+scale_fill_manual(values =c("#555555","#c699c6"))+scale_color_manual(values =c("#555555","#c699c6"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
p2 + facet_wrap( ~ Patient, scales="free",ncol=4)



############# PERFORM CELLCHAT V2 FOR RECEPTOR LIGAND INTERACTION ANALYSIS (https://github.com/jinworks/CellChat) ################
#follow this vignette: https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

BM.combined.sct <- LoadSeuratRds("BM.combined.sct_2025.rds")

library(DescTools)

BM.combined.sct$SampleInfo_Type_Time <- paste(BM.combined.sct$Patient,BM.combined.sct$SampleType,BM.combined.sct$Time,sep="-")
BM.combined.sct$SampleInfo_Type_Time<- factor(BM.combined.sct$SampleInfo_Type_Time, levels = c("24780-BM-Scr","24780-BM-PC2","122759-BM-Scr","122759-BM-PC2","113446-BM-Scr","113446-BM-PC2","122739-BM-Scr","122739-BM-PC2","121535-BM-Scr","121535-BM-PC2","121474-BM-Scr","121474-BM-PC2","119837-BM-Scr","119837-BM-PC2","121015-BM-Scr","121015-BM-PC2","122976-BM-Scr","122976-BM-PC2","114152-BM-Scr","114152-BM-PC2","119325-BM-Scr","119325-BM-PC2","118877-BM-Scr","118877-BM-PC2","115585-BM-Scr","115585-BM-PC2","24780-PBMC-Scr","24780-PBMC-C2D1","24780-PBMC-C2D8","121474-PBMC-Scr","121474-PBMC-C1D8","121474-PBMC-C2D1","121015-PBMC-Scr","121015-PBMC-C1D8","121015-PBMC-C2D1","118877-PBMC-Scr","118877-PBMC-C1D8","118877-PBMC-C2D1"))

BM.combined.sct@meta.data$PatientTime <- paste(BM.combined.sct@meta.data$Patient,BM.combined.sct@meta.data$Time2, sep="_")
BM.combined.sct@meta.data$TimeSurvival <- paste(BM.combined.sct@meta.data$Survival,BM.combined.sct@meta.data$Time2,sep = "_")

BM.combined.sct$SampleInfo_Type_Time_Specific <- paste(BM.combined.sct$manualAnno_specific2,BM.combined.sct$SampleInfo_Type_Time, sep ="_")
BM.combined.sct$SampleInfo_Type_Time_broad <- paste(BM.combined.sct$manualAnno_specific,BM.combined.sct$SampleInfo_Type_Time, sep ="_")

BM.combined.sct.BM <- subset(BM.combined.sct, subset = SampleType == "BM")

Idents(BM.combined.sct.BM) <- "manualAnno_specific2"  

samplelistname <- data.frame()
df.net_ALL <-data.frame()

for (i in c("Short_Pre","Short_Post","Long_Pre","Long_Post")) {
print(i)
###Prepare required input data for CellChat analysis per patient
bm <- subset(BM.combined.sct.BM, idents = c("Stromal","MAIT", "Megakaryocyte","Differentiating Stem Cell","Neutrophil"), invert = TRUE) #remove low cell count immune types

Idents(bm) <- "manualAnno_specific2"  
bm1 <- subset(bm, subset = TimeSurvival == i)

for (j in c(unique(BM.combined.sct.BM@meta.data[BM.combined.sct.BM@meta.data$TimeSurvival == i,]$Patient))) {
print(j)
bm2 <- subset(bm1, subset = Patient == j) #perform patient specific 

input <- GetAssayData(object = bm2, assay = 'SCT', slot = "data")
labels <- Idents(bm2)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe


bm2 <- NULL

###Create a CellChat object
cellChat <- createCellChat(object = input, meta = meta, group.by = "labels")

###Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)

# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 

# set the used database in the object
cellChat@DB <- CellChatDB.use

tryCatch({
###Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellChat <- subsetData(cellChat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellChat <- identifyOverExpressedGenes(cellChat)
#An object of class CellChat created from a single dataset 
#22137 genes.
#61487 cells. 
#CellChat analysis of single cell RNA-seq data! 

cellChat <- identifyOverExpressedInteractions(cellChat)
#The number of highly variable ligand-receptor pairs used for signaling inference is 699 

###Compute the communication probability and infer cellular communication network
cellChat <- computeCommunProb(cellChat, type = "triMean")

cellChat <- filterCommunication(cellChat, min.cells = 10) #similar to pseudobulk cutoff

#Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellChat) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
df.net <- df.net[df.net$pval <= 0.5 & df.net$pval > 0,]
df.net$Comparison <- i
df.net$Patient <- j
df.net_ALL <- rbind(df.net_ALL, df.net)

df.net <- NULL

#Infer the cell-cell communication at a signaling pathway level
cellChat <- computeCommunProbPathway(cellChat)


#Calculate the aggregated cell-cell communication network
cellChat <- aggregateNet(cellChat)
assign(paste("cellChat",i,j, sep = "_"), cellChat)

sampleinfo <- data.frame("Name" = c(paste("cellChat",i,j, sep = "_")), "Comparison" = c(i), "Patient" = c(j)) 
print(sampleinfo)
samplelistname <- rbind(samplelistname,sampleinfo)
  
cellChat <- NULL},
error=function(e){cat("ERROR :",conditionMessage(e), "\n")})}}



#combine all cells into SurvivalTime
options(future.globals.maxSize= 1248576000)

df.net_ALL_total <-data.frame()

for (i in c("Short_Pre","Short_Post","Long_Pre","Long_Post")) {
  print(i)
  ###Prepare required input data for CellChat analysis per patient
  bm <- subset(BM.combined.sct.BM, idents = c("Stromal","MAIT", "Megakaryocyte","Differentiating Stem Cell","Neutrophil"), invert = TRUE) #remove low cell count immune types
  
  Idents(bm) <- "manualAnno_specific2"  
  bm <- subset(bm, subset = TimeSurvival == i)

    input <- GetAssayData(object = bm, assay = 'SCT', slot = "data")
    labels <- Idents(bm)
    meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe
    
    bm <- NULL
    
    ###Create a CellChat object
    cellChat <- createCellChat(object = input, meta = meta, group.by = "labels")
    
    ###Set the ligand-receptor interaction database
    CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
    showDatabaseCategory(CellChatDB)

    # use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
    CellChatDB.use <- subsetDB(CellChatDB)
    
    # set the used database in the object
    cellChat@DB <- CellChatDB.use
    
    tryCatch({
    ###Preprocessing the expression data for cell-cell communication analysis
    # subset the expression data of signaling genes for saving computation cost
    cellChat <- subsetData(cellChat) # This step is necessary even if using the whole database
    future::plan("multisession", workers = 4) # do parallel
    cellChat <- identifyOverExpressedGenes(cellChat)
    #An object of class CellChat created from a single dataset 
    #22137 genes.
    #61487 cells. 
    #CellChat analysis of single cell RNA-seq data! 
    
    cellChat <- identifyOverExpressedInteractions(cellChat)
    #The number of highly variable ligand-receptor pairs used for signaling inference is 699 
    
    # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
    #cellChat2 <- projectData(cellChat, PPI.human)
    
    
    ###Compute the communication probability and infer cellular communication network
    cellChat <- computeCommunProb(cellChat, type = "triMean")
    #cellChat <- computeCommunProb(cellChat, type =  "truncatedMean", trim = 0.1) 
    
    cellChat <- filterCommunication(cellChat, min.cells = 20) #similar to pseudobulk cutoff
    
    #Extract the inferred cellular communication network as a data frame
    df.net <- subsetCommunication(cellChat) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
    df.net <- df.net[df.net$pval <= 0.5 & df.net$pval > 0,]
    df.net$Comparison <- i
    df.net_ALL_total <- rbind(df.net_ALL_total, df.net)
    
    df.net <- NULL
    
    #df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
    
    #df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) #gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
    
    #Infer the cell-cell communication at a signaling pathway level
    cellChat <- computeCommunProbPathway(cellChat)
    
    
    #Calculate the aggregated cell-cell communication network
    cellChat <- aggregateNet(cellChat)
    assign(paste("cellChat",i, sep = "_"), cellChat)

    cellChat <- NULL},
    error=function(e){cat("ERROR :",conditionMessage(e), "\n")})}




#make sure cells are in right order
cellChat_Short_Pre <- updateClusterLabels(cellChat_Short_Pre, new.order = c("EPC","MEP","HSC","LMPP","GMP","Classical Mono","Intermediate Mono","Non-classical Mono","mDC","pDC","Naive B","Memory B","Plasmablast","Vd1 gdT","Vd2 gdT","CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL","CD4 Treg","CD8 Naive","CD8 TCM","CD8 TEM","CD8 TEMRA","CD56dim NK","CD56bright NK"))
cellChat_Short_Post <- updateClusterLabels(cellChat_Short_Post, new.order = c("EPC","MEP","HSC","LMPP","GMP","Classical Mono","Intermediate Mono","Non-classical Mono","mDC","pDC","Naive B","Memory B","Plasmablast","Vd1 gdT","Vd2 gdT","CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL","CD4 Treg","CD8 Naive","CD8 TCM","CD8 TEM","CD8 TEMRA","CD56dim NK","CD56bright NK"))
cellChat_Long_Pre <- updateClusterLabels(cellChat_Long_Pre, new.order = c("EPC","MEP","HSC","LMPP","GMP","Classical Mono","Intermediate Mono","Non-classical Mono","mDC","pDC","Naive B","Memory B","Plasmablast","Vd1 gdT","Vd2 gdT","CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL","CD4 Treg","CD8 Naive","CD8 TCM","CD8 TEM","CD8 TEMRA","CD56dim NK","CD56bright NK"))
cellChat_Long_Post <- updateClusterLabels(cellChat_Long_Post, new.order = c("EPC","MEP","HSC","LMPP","GMP","Classical Mono","Intermediate Mono","Non-classical Mono","mDC","pDC","Naive B","Memory B","Plasmablast","Vd1 gdT","Vd2 gdT","CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL","CD4 Treg","CD8 Naive","CD8 TCM","CD8 TEM","CD8 TEMRA","CD56dim NK","CD56bright NK"))


#######Comparison analysis of multiple datasets using CellChat ###########
object.list <- list(ShortPre = cellChat_Short_Pre, ShortPost = cellChat_Short_Post)
cellchat_ShortPrePost <- mergeCellChat(object.list, add.names = names(object.list))
object.list <- list(LongPre = cellChat_Long_Pre, LongPost = cellChat_Long_Post)
cellchat_LongPrePost <- mergeCellChat(object.list, add.names = names(object.list))
object.list <- list(ShortPre = cellChat_Short_Pre, LongPre = cellChat_Long_Pre)
cellchat_ShortLongPre <- mergeCellChat(object.list, add.names = names(object.list))
object.list <- list(ShortPost = cellChat_Short_Post, LongPost = cellChat_Long_Post)
cellchat_ShortLongPost <- mergeCellChat(object.list, add.names = names(object.list))



#Short Pre vs. Post
gg1 <- compareInteractions(cellchat_ShortPrePost, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_ShortPrePost, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#Long Pre vs. Post
gg1 <- compareInteractions(cellchat_LongPrePost, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_LongPrePost, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#Short Pre vs. Long Pre
gg1 <- compareInteractions(cellchat_ShortLongPre, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_ShortLongPre, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#Short Post vs. Long Post
gg1 <- compareInteractions(cellchat_ShortLongPost, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_ShortLongPost, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2



########Start analaysis of Short Pre vs Short Post #####

#(A) Circle plot showing differential number of interactions or interaction strength among different cell populations across two datasets
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat_ShortPrePost, weight.scale = T)
netVisual_diffInteraction(cellchat_ShortPrePost, weight.scale = T, measure = "weight")


#(B) Heatmap showing differential number of interactions or interaction strength among different cell populations across two datasets
gg1 <- netVisual_heatmap(cellchat_ShortPrePost)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat_ShortPrePost, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2


##Compare the major sources and targets in a 2D space
#(A) Identify cell populations with significant changes in sending or receiving signals
object.list <- list(ShortPre = cellChat_Short_Pre, ShortPost = cellChat_Short_Post,LongPre = cellChat_Long_Pre, LongPost = cellChat_Long_Post) #all four datasets
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)


#(B) Identify the signaling changes of specific cell populations
gg1 <-netAnalysis_signalingChanges_scatter(cellchat, idents.use = "HSC", comparison = c(1,3), signaling.exclude = c("MIF","GALECTIN","SELPLG","MHC-I","MHC-II","CD99","CypA","MK"))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "LMPP", comparison = c(1,3), signaling.exclude = c("MIF","GALECTIN","SELPLG","MHC-I","MHC-II","CD99","CypA","MK"))
patchwork::wrap_plots(plots = list(gg1,gg2))

gg1 <-netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD4 Naive", comparison = c(1,3), signaling.exclude = c("MIF","GALECTIN","SELPLG","MHC-I","MHC-II","CD99","CypA","MK"))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD4 TCM", comparison = c(1,3), signaling.exclude = c("MIF","GALECTIN","SELPLG","MHC-I","MHC-II","CD99","CypA","MK"))
patchwork::wrap_plots(plots = list(gg1,gg2))



#Part II: Identify altered signaling with distinct network architecture and interaction strength
#Identify signaling networks with larger (or less) difference as well as signaling groups based on their functional/structure similarity
#Functional similarity: High degree of functional similarity indicates major senders and receivers are similar, and it can be interpreted as the two signaling pathways or two ligand-receptor pairs exhibit similar and/or redundant roles. NB: Functional similarity analysis is not applicable to multiple datsets with different cell type composition.
#Structural similarity: A structural similarity was used to compare their signaling network structure, without considering the similarity of senders and receivers. NB: Structural similarity analysis is applicable to multiple datsets with the same cell type composition or the vastly different cell type composition.
#since our samples dont have similar cell composition, we will use structural similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural",comparison = c(1,2))
cellchat <- netEmbedding(cellchat, type = "structural",comparison = c(1,2))
cellchat <- netClustering(cellchat, type = "structural",comparison = c(1,2))
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5,comparison = c(1,2))

cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural",comparison = c(1,3))
cellchat <- netEmbedding(cellchat, type = "structural",comparison = c(1,3))
cellchat <- netClustering(cellchat, type = "structural",comparison = c(1,3))
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5,comparison = c(1,3))

cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural",comparison = c(3,4))
cellchat <- netEmbedding(cellchat, type = "structural",comparison = c(3,4))
cellchat <- netClustering(cellchat, type = "structural",comparison = c(3,4))
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5,comparison = c(3,4))

cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural",comparison = c(2,4))
cellchat <- netEmbedding(cellchat, type = "structural",comparison = c(2,4))
cellchat <- netClustering(cellchat, type = "structural",comparison = c(2,4))
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5,comparison = c(2,4))




#Identify altered signaling with distinct interaction strength
#(A) Compare the overall information flow of each signaling pathway or ligand-receptor pair

gg6 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE,comparison = c(1,3)) #short pre vs long pre
gg8 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE,comparison = c(2,4)) #short post vs short long

gg6 + gg8

#(B) Compare outgoing (or incoming) signaling patterns associated with each cell population
library(ComplexHeatmap)

# combining all the identified signaling pathways from different datasets 

# Short Pre vs Long Pre
pathway.union <- union(object.list[[1]]@netP$pathways, object.list[[3]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "all", signaling = pathway.union, title = names(object.list)[1], width = 9, height = 18)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[3]], pattern = "all", signaling = pathway.union, title = names(object.list)[3], width = 9, height = 18)

draw(ht1 + ht2, ht_gap = unit(1, "cm"))

# Short Post vs Long Post
pathway.union <- union(object.list[[2]]@netP$pathways, object.list[[4]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "all", signaling = pathway.union, title = names(object.list)[2], width = 9, height = 18)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[4]], pattern = "all", signaling = pathway.union, title = names(object.list)[4], width = 9, height = 18)

draw(ht1 + ht2, ht_gap = unit(1, "cm"))





############### PSEUDO BULK ANALYSIS ON PBMC SAMPLES ########################
#START NEW SESSION TO RESET MEMORY CONSUMPTION FROM previous DEG calling
BM.combined.sct <- LoadSeuratRds("BM.combined.sct_2025.rds") ######## LOAD ONLY THE BM SEURAT DATASET TO REDUCE MEMORY CONSUMPTION AND SERVER CRASH

BM.combined.sct$SampleInfo_Type_Time <- paste(BM.combined.sct$Patient,BM.combined.sct$SampleType,BM.combined.sct$Time,sep="-")
BM.combined.sct$SampleInfo_Type_Time<- factor(BM.combined.sct$SampleInfo_Type_Time, levels = c("24780-BM-Scr","24780-BM-PC2","122759-BM-Scr","122759-BM-PC2","113446-BM-Scr","113446-BM-PC2","122739-BM-Scr","122739-BM-PC2","121535-BM-Scr","121535-BM-PC2","121474-BM-Scr","121474-BM-PC2","119837-BM-Scr","119837-BM-PC2","121015-BM-Scr","121015-BM-PC2","122976-BM-Scr","122976-BM-PC2","114152-BM-Scr","114152-BM-PC2","119325-BM-Scr","119325-BM-PC2","118877-BM-Scr","118877-BM-PC2","115585-BM-Scr","115585-BM-PC2","24780-PBMC-Scr","24780-PBMC-C2D1","24780-PBMC-C2D8","121474-PBMC-Scr","121474-PBMC-C1D8","121474-PBMC-C2D1","121015-PBMC-Scr","121015-PBMC-C1D8","121015-PBMC-C2D1","118877-PBMC-Scr","118877-PBMC-C1D8","118877-PBMC-C2D1"))

BM.combined.sct@meta.data$PatientTime <- paste(BM.combined.sct@meta.data$Patient,BM.combined.sct@meta.data$Time2, sep="_")
BM.combined.sct@meta.data$TimeSurvival <- paste(BM.combined.sct@meta.data$Survival,BM.combined.sct@meta.data$Time2,sep = "_")

BM.combined.sct.BM <- subset(BM.combined.sct, subset = SampleType == "BM")
BM.combined.sct.PBMC <- subset(BM.combined.sct, subset = SampleType == "PBMC")
DefaultAssay(BM.combined.sct.PBMC) = "RNA"
BM.combined.sct.PBMC <- JoinLayers(BM.combined.sct.PBMC)

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- BM.combined.sct.PBMC[["RNA"]]$counts 

# Set up metadata as desired for aggregation and DE analysis
BM.combined.sct.PBMC$sample_id <- paste(BM.combined.sct.PBMC$Patient,BM.combined.sct.PBMC$TimeSurvival, sep = "_")

BM.combined.sct.PBMC$sample_id <- factor(BM.combined.sct.PBMC$sample_id)


metadata <- BM.combined.sct.PBMC@meta.data


# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)


# Single-cell RNA-seq analysis - Pseudobulk DE analysis with DESeq2library(tidyverse)
library(cowplot)
library(edgeR)
library(Matrix)
library(Matrix.utils) #remotes::install_github("cvarrichio/Matrix.utils")
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)

# Extract unique names of clusters (= levels of cluster_id factor variable)
cluster_names_broad <- unique(colData(sce)$manualAnno_specific)
cluster_names_broad 

cluster_names_specific <- unique(colData(sce)$manualAnno_specific2)
cluster_names_specific

# Extract unique names of samples (= levels of sample_id factor variable)
sample_names <- levels(colData(sce)$sample_id)
sample_names


# Subset metadata to include only the variables you want to aggregate across (here, we want to aggregate by sample and by cluster)
groups_broad <- colData(sce)[, c("manualAnno_specific","sample_id")]
head(groups_broad)

groups_specific <- colData(sce)[, c("manualAnno_specific2","sample_id")]
head(groups_specific)


# Aggregate across cluster-sample groups
# transposing row/columns to have cell_ids as row names matching those of groups
aggr_counts_broad <- aggregate.Matrix(t(counts(sce)), 
                                      groupings = groups_broad, fun = "sum") 
aggr_counts_specific <- aggregate.Matrix(t(counts(sce)), 
                                         groupings = groups_specific, fun = "sum") 


# Transpose aggregated matrix to have genes as rows and samples as columns
aggr_counts_broad <- t(aggr_counts_broad)
aggr_counts_specific <- t(aggr_counts_specific)



# Loop over all cell types to extract corresponding counts, and store information in a list

## Initiate empty list
counts_ls_broad <- list()

for (i in 1:length(cluster_names_broad)) {
  
  ## Extract indexes of columns in the global matrix that match a given cluster
  column_idx <- which(tstrsplit(colnames(aggr_counts_broad), "_")[[1]] == cluster_names_broad[i])
  
  ## Store corresponding sub-matrix as one element of a list
  counts_ls_broad[[i]] <- aggr_counts_broad[, column_idx]
  names(counts_ls_broad)[i] <- as.character(cluster_names_broad[i])
  
}

# Explore the different components of the list
str(counts_ls_broad)



## Initiate empty list
counts_ls_specific <- list()

for (i in 1:length(cluster_names_specific)) {
  
  ## Extract indexes of columns in the global matrix that match a given cluster
  column_idx <- which(tstrsplit(colnames(aggr_counts_specific), "_")[[1]] == cluster_names_specific[i])
  
  ## Store corresponding sub-matrix as one element of a list
  counts_ls_specific[[i]] <- aggr_counts_specific[, column_idx]
  names(counts_ls_specific)[i] <- as.character(cluster_names_specific[i])
  
}

# Explore the different components of the list
str(counts_ls_specific)


# Reminder: explore structure of metadata
head(colData(sce))

# Extract sample-level variables
metadata <- colData(sce) %>% 
  as.data.frame() %>% 
  dplyr::select(TimeSurvival,Patient,sample_id)

dim(metadata)
head(metadata)

# Exclude duplicated rows
metadata <- metadata[!duplicated(metadata), ]

dim(metadata)
head(metadata)

# Rename rows
rownames(metadata) <- metadata$sample_id
head(metadata)


# Number of cells per sample and cluster
t <- table(colData(sce)$sample_id,colData(sce)$manualAnno_specific)
t2 <- table(colData(sce)$sample_id,colData(sce)$manualAnno_specific2)

t ### 
B cell CD34  CD4  CD8   DC Differentiating Stem Cell  EPC  gdT MAIT Megakaryocyte Monocyte Neutrophil   NK Plasmablast
118877_Long_C1D8    190    1 1482  414    1                        11    0   45   58             6       22          1  517           0
118877_Long_Post    269   16 1375  162    3                        60    8   10    6            16      192          2   43           4
118877_Long_Pre     195    3 1577  338    1                         8    0   38   50            12       15          0  392           4
121015_Long_C1D8     40    4  495  364    0                        21    2   48    1             4      118        119   53           0
121015_Long_Post     39   42  748  624   17                        41    1   59    3            21      117        137  116           0
121015_Long_Pre      61   13  412  185   12                        22    0   17    1             5      130         82   26           4
121474_Long_C1D8     31   42  359   64    3                        20    0   21   16           135       46          0  260           1
121474_Long_Post     39   48  561   84    0                        20    7    9   16            12        2          1   54          21
121474_Long_Pre      84  157 1071  117    8                        47    1   42   92             7       16          0  676           2
24780_Short_C2D8    286   31 2026  955    1                        44    2  444  109             5        0          0  430           3
24780_Short_Post    293   62  434  225    0                        46    0   78   18            18        1          0   47          36
24780_Short_Pre     135   41 1440  422    2                         5    0  175   94             2        4          0  190           4

t2 ##                    CD4 CTL CD4 Naive CD4 TCM CD4 TEM CD4 Treg CD56bright NK CD56dim NK CD8 Naive CD8 TCM CD8 TEM CD8 TEMRA Classical Mono Differentiating Stem Cell  EPC  GMP  HSC Intermediate Mono LMPP MAIT  mDC Megakaryocyte Memory B  MEP Naive B Neutrophil Non-classical Mono  pDC Plasmablast Vd1 gdT Vd2 gdT
118877_Long_C1D8     115       736     606       3       22            19        498        40      42     107       225             13                        11    0    0    1                 6    0   58    0             6      146    0      44          1                  3    1           0      45       0
118877_Long_Post      27       771     551       1       25             2         41        44      19      42        57            144                        60    8    0    5                33    4    6    2            16      210    7      59          2                 15    1           4       9       1
118877_Long_Pre       66       818     665       0       28            24        368        49      35      96       158              9                         8    0    0    2                 1    1   50    1            12      146    0      49          0                  5    0           4      33       5
121015_Long_C1D8      12        87     385       0       11             0         53         4       8      68       284            101                        21    2    1    0                15    3    1    0             4       31    0       9        119                  2    0           0      44       4
121015_Long_Post      23       301     401       0       23             7        109         7       8     250       359             63                        41    1    0   24                44   13    3   16            21       25    5      14        137                 10    1           0      55       4
121015_Long_Pre        3       103     291       1       14             2         24         6       9      39       131             97                        22    0    1    6                22    6    1   11             5       49    1      12         82                 10    1           4      15       2
121474_Long_C1D8       4       125     221       3        6            13        247         8       4      49         3             19                        20    0    1   22                15   13   16    2           135       24    6       7          0                 12    1           1       0      21
121474_Long_Post       2       225     293       0       41            13         41        15       7      39        23              2                        20    7    0   20                 0    8   16    0            12       17   20      22          1                  0    0          21       2       7
121474_Long_Pre        2       324     722       5       18           208        468        33      20      48        16              5                        47    1    2   96                 2   33   92    2             7       51   26      33          0                  9    6           2       4      38
24780_Short_C2D8      82       466    1462       1       15            11        419        64      47     307       537              0                        44    2    7    8                 0   15  109    0             5      274    1      12          0                  0    1           3     142     302
24780_Short_Post      38        31     346       0       19             4         43         3      11     106       105              0                        46    0   24   23                 0   13   18    0            18      282    2      11          0                  1    0          36      10      68
24780_Short_Pre      101       465     844       0       30            10        180        25      32     169       196              2                         5    0    7   22                 0   10   94    0             2      126    2       9          0                  2    2           4      41     134

###NOTE THAT LOT OF CELLS ARE LOW ABUNDANCE PER PATIENT PER TIMEPOINT. Enough cells for CD4 Naive, CD4 TCM, CD56dim NK, CD8 Naive, CD8 TEM, CD8 TEMRA, Classical Mono, GMP, HSC, LMPP,Vd1 gdT


# Creating metadata list
## Initiate empty list
metadata_ls_broad_PBMC <- list()

for (i in c("CD34","Monocyte","B cell","CD4","CD8","NK","gdT")) {
  
  ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
  df <- data.frame(cluster_sample_id = colnames(counts_ls_broad[[i]]))
  
  ## Use tstrsplit() to separate cluster (cell type) and sample IDs
  df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
  df$sample_id  <- paste(tstrsplit(df$cluster_sample_id, "_")[[2]],tstrsplit(df$cluster_sample_id, "_")[[3]],tstrsplit(df$cluster_sample_id, "_")[[4]],sep="_")
  
  
  ## Retrieve cell count information for this cluster from global cell count table
  idx <- which(colnames(t) == unique(df$cluster_id))
  cell_counts <- t[, idx]
  
  ## Remove samples with less than 10 cell contributing to the cluster
  cell_counts <- cell_counts[cell_counts > 10]
  
  ## Match order of cell_counts and sample_ids
  sample_order <- match(df$sample_id, names(cell_counts))
  cell_counts <- cell_counts[sample_order]
  
  ## Append cell_counts to data frame
  df$cell_count <- cell_counts
  
  
  ## Join data frame (capturing metadata specific to cluster) to generic metadata
  df <- plyr::join(df, metadata, 
                   by = intersect(names(df), names(metadata)))
  
  ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
  rownames(df) <- df$cluster_sample_id
  
  #remove samples that don't pass cell # cutoff
  df <- df[complete.cases(df),]
  
  ## Store complete metadata for cluster i in list
  metadata_ls_broad_PBMC[[i]] <- df
}

# Explore the different components of the list
str(metadata_ls_broad_PBMC)




## Initiate empty list
metadata_ls_specific_PBMC <- list()

for (i in 1:length(counts_ls_specific)) {
  
  ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
  df <- data.frame(cluster_sample_id = colnames(counts_ls_specific[[i]]))
  
  ## Use tstrsplit() to separate cluster (cell type) and sample IDs
  df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
  df$sample_id  <- paste(tstrsplit(df$cluster_sample_id, "_")[[2]],tstrsplit(df$cluster_sample_id, "_")[[3]],tstrsplit(df$cluster_sample_id, "_")[[4]],sep="_")
  
  
  ## Retrieve cell count information for this cluster from global cell count table
  idx <- which(colnames(t2) == unique(df$cluster_id))
  cell_counts <- t2[, idx]
  
  ## Remove samples with < 10 cell contributing to the cluster
  cell_counts <- cell_counts[cell_counts > 10]
  
  ## Match order of cell_counts and sample_ids
  sample_order <- match(df$sample_id, names(cell_counts))
  cell_counts <- cell_counts[sample_order]
  
  ## Append cell_counts to data frame
  df$cell_count <- cell_counts
  
  
  ## Join data frame (capturing metadata specific to cluster) to generic metadata
  df <- plyr::join(df, metadata, 
                   by = intersect(names(df), names(metadata)))
  
  ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
  rownames(df) <- df$cluster_sample_id
  
  ## Store complete metadata for cluster i in list
  metadata_ls_specific_PBMC[[i]] <- df
  names(metadata_ls_specific_PBMC)[i] <- unique(df$cluster_id)
  metadata_ls_specific_PBMC[[i]] <- metadata_ls_specific_PBMC[[i]][complete.cases(metadata_ls_specific_PBMC[[i]]),]
  print(unique(df$cluster_id))
  print(table(metadata_ls_specific_PBMC[[i]]$TimeSurvival))
}




##generate normalized counts file
DESeq2_NormCounts_broad_PBMC <-data.frame()
for (broadcelltype in c("CD34","Monocyte","B cell","CD4","CD8","NK","gdT")){
  print(broadcelltype)
  cluster_counts <- counts_ls_broad[[broadcelltype]]
  cluster_metadata <- metadata_ls_broad_PBMC[[broadcelltype]]
  cluster_metadata$Time <- apply(cluster_metadata, 1, function(x) {unlist(strsplit(x[5],"_"))[2]})
  cluster_metadata$Patient <- factor(cluster_metadata$Patient)
  
  cluster_counts <- cluster_counts[,cluster_metadata$cluster_sample_id]
  
  print(all(colnames(cluster_counts) == rownames(cluster_metadata)))
  
  # Create DESeq2 object    
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ Time)
  keep <- rowSums(counts(dds)) >= 1 #filter non expressed genes (not doing GSEA)
  dds <- dds[keep,]#filter lowly expressed genes
  dds <- DESeq(dds)
  
  NormCountsData <- data.frame(counts(dds, normalize = T), stringsAsFactors = F)
  NormCountsData$Gene <- rownames(NormCountsData)
  if (broadcelltype == "CD34"){
    d <- NormCountsData
    DESeq2_NormCounts_broad_PBMC <- NormCountsData }
  else {
    d2 <- NormCountsData
    DESeq2_NormCounts_broad_PBMC <- merge(DESeq2_NormCounts_broad_PBMC,NormCountsData,by = "Gene", all = T) }
}





DESeq2_NormCounts_specific_PBMC <-data.frame()
for (specificcelltype in c("CD4 Naive", "CD4 TCM", "CD56dim NK", "CD8 Naive", "CD8 TEM", "CD8 TEMRA", "Naive B","Memory B","Classical Mono", "HSC", "LMPP","Vd1 gdT","Vd2 gdT")){
  tryCatch({ print(specificcelltype)
    cluster_counts <- counts_ls_specific[[specificcelltype]]
    cluster_metadata <- metadata_ls_specific_PBMC[[specificcelltype]]
    cluster_metadata$Time <- apply(cluster_metadata, 1, function(x) {unlist(strsplit(x[5],"_"))[2]})
    cluster_metadata$Patient <- factor(cluster_metadata$Patient)
    
    cluster_counts <- cluster_counts[,cluster_metadata$cluster_sample_id]
    
    print(all(colnames(cluster_counts) == rownames(cluster_metadata)))
    
    
    # Create DESeq2 object    
    dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                  colData = cluster_metadata, 
                                  design = ~ Time)
    keep <- rowSums(counts(dds)) >= 10 #DONT FILTER TO CAPTURE ALL GENES (not doing GSEA)
    dds <- dds[keep,]#filter lowly expressed genes
    dds <- DESeq(dds)
    
    NormCountsData <- data.frame(counts(dds, normalize = T), stringsAsFactors = F)
    NormCountsData$Gene <- rownames(NormCountsData)
    if (specificcelltype == "CD4 Naive"){
      d <- NormCountsData
      DESeq2_NormCounts_specific_PBMC <- NormCountsData }
    else {
      d2 <- NormCountsData
      DESeq2_NormCounts_specific_PBMC <- merge(DESeq2_NormCounts_specific_PBMC,NormCountsData,by = "Gene", all = T) }
  })}

#save.image("BM_CITEseq_DESeq2_PseudobulkPBMC_seurat5_051424.RData")




### plot pseudobulk 

GOIgenes <- c("CDKN1A","CDKN1B","CDKN2C","CDKN2D","CXCL8","CXCL2","CXCL3") #CD34

##broad anno
DESeq2_NormCounts_broad_PBMC_GOI <- DESeq2_NormCounts_broad_PBMC[DESeq2_NormCounts_broad_PBMC$Gene %in% GOIgenes,]

DESeq2_NormCounts_broad_PBMC_GOI1 <- DESeq2_NormCounts_broad_PBMC_GOI
DESeq2_NormCounts_broad_PBMC_GOI1[is.na(DESeq2_NormCounts_broad_PBMC_GOI1)] <- -0.1
DESeq2_NormCounts_broad_PBMC_GOI2 <- DESeq2_NormCounts_broad_PBMC_GOI1[,c("Gene","CD34_24780_Short_Pre","CD34_24780_Short_Post","CD34_24780_Short_C2D8", "CD34_121474_Long_Pre","CD34_121474_Long_C1D8","CD34_121474_Long_Post","CD34_121015_Long_Pre","CD34_121015_Long_Post","CD34_118877_Long_Post","Monocyte_121474_Long_Pre","Monocyte_121474_Long_C1D8","Monocyte_121015_Long_Pre","Monocyte_121015_Long_C1D8","Monocyte_121015_Long_Post","Monocyte_118877_Long_Pre","Monocyte_118877_Long_C1D8","Monocyte_118877_Long_Post","B.cell_24780_Short_Pre","B.cell_24780_Short_Post","B.cell_24780_Short_C2D8", "B.cell_121474_Long_Pre","B.cell_121474_Long_C1D8","B.cell_121474_Long_Post", "B.cell_121015_Long_Pre","B.cell_121015_Long_C1D8","B.cell_121015_Long_Post", "B.cell_118877_Long_Pre","B.cell_118877_Long_C1D8","B.cell_118877_Long_Post","CD4_24780_Short_Pre","CD4_24780_Short_Post","CD4_24780_Short_C2D8", "CD4_121474_Long_Pre","CD4_121474_Long_C1D8","CD4_121474_Long_Post", "CD4_121015_Long_Pre","CD4_121015_Long_C1D8","CD4_121015_Long_Post", "CD4_118877_Long_Pre","CD4_118877_Long_C1D8","CD4_118877_Long_Post","CD8_24780_Short_Pre","CD8_24780_Short_Post","CD8_24780_Short_C2D8", "CD8_121474_Long_Pre","CD8_121474_Long_C1D8","CD8_121474_Long_Post", "CD8_121015_Long_Pre","CD8_121015_Long_C1D8","CD8_121015_Long_Post", "CD8_118877_Long_Pre","CD8_118877_Long_C1D8","CD8_118877_Long_Post","NK_24780_Short_Pre","NK_24780_Short_Post","NK_24780_Short_C2D8", "NK_121474_Long_Pre","NK_121474_Long_C1D8","NK_121474_Long_Post", "NK_121015_Long_Pre","NK_121015_Long_C1D8","NK_121015_Long_Post", "NK_118877_Long_Pre","NK_118877_Long_C1D8","NK_118877_Long_Post","gdT_24780_Short_Pre","gdT_24780_Short_Post","gdT_24780_Short_C2D8", "gdT_121474_Long_Pre","gdT_121474_Long_C1D8", "gdT_121015_Long_Pre","gdT_121015_Long_C1D8","gdT_121015_Long_Post", "gdT_118877_Long_Pre","gdT_118877_Long_C1D8")]



#look at CD34 specifically
DESeq2_NormCounts_broad_PBMC_GOI2 <- DESeq2_NormCounts_broad_PBMC_GOI2[,c("Gene","CD34_24780_Short_Pre","CD34_24780_Short_Post","CD34_24780_Short_C2D8", "CD34_121474_Long_Pre","CD34_121474_Long_C1D8","CD34_121474_Long_Post","CD34_121015_Long_Pre","CD34_121015_Long_Post","CD34_118877_Long_Post")]

library("circlize") ## For color options
name <- DESeq2_NormCounts_broad_PBMC_GOI2[,c("Gene")]
df.OG2 <- data.matrix(log2(DESeq2_NormCounts_broad_PBMC_GOI2[,c(2:10)]+0.01))
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(DESeq2_NormCounts_broad_PBMC_GOI2[,c(2:10)])
pheatmap(df.OG2,color = colorRampPalette(c("#f0f0f0",rev(brewer.pal(n = 11, name = "Spectral"))))(100),breaks = c(seq(-0.1,15,by=0.151)),show_rownames = T,cluster_rows =F, cluster_cols=F)

## Generate Z-score for TPM from all samples (USE THIS)###
library(matrixStats)
DESeq2_NormCounts_broad_PBMC_GOI21<-log2(DESeq2_NormCounts_broad_PBMC_GOI2[,c(2:10)]+0.01)
DESeq2_NormCounts_broad_PBMC_GOI21_Zscore<- (DESeq2_NormCounts_broad_PBMC_GOI21-rowMeans(DESeq2_NormCounts_broad_PBMC_GOI21))/(rowSds(as.matrix(DESeq2_NormCounts_broad_PBMC_GOI21)))[row(DESeq2_NormCounts_broad_PBMC_GOI21)]
name <- DESeq2_NormCounts_broad_PBMC_GOI2[,c("Gene")]
df.OG2 <- data.matrix(DESeq2_NormCounts_broad_PBMC_GOI21_Zscore)
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(DESeq2_NormCounts_broad_PBMC_GOI2[,c(2:10)])
pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = T,cluster_rows =F, cluster_cols=F)
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),clustering_method="ward.D2",show_rownames = T,cluster_rows =T, cluster_cols=F)





GOIgenes <- c("S100A8","S100A9") #Monocyte

DESeq2_NormCounts_broad_PBMC_GOI <- DESeq2_NormCounts_broad_PBMC[DESeq2_NormCounts_broad_PBMC$Gene %in% GOIgenes,]

DESeq2_NormCounts_broad_PBMC_GOI1 <- DESeq2_NormCounts_broad_PBMC_GOI
DESeq2_NormCounts_broad_PBMC_GOI1[is.na(DESeq2_NormCounts_broad_PBMC_GOI1)] <- -0.1
DESeq2_NormCounts_broad_PBMC_GOI2 <- DESeq2_NormCounts_broad_PBMC_GOI1[,c("Gene","CD34_24780_Short_Pre","CD34_24780_Short_Post","CD34_24780_Short_C2D8", "CD34_121474_Long_Pre","CD34_121474_Long_C1D8","CD34_121474_Long_Post","CD34_121015_Long_Pre","CD34_121015_Long_Post","CD34_118877_Long_Post","Monocyte_121474_Long_Pre","Monocyte_121474_Long_C1D8","Monocyte_121015_Long_Pre","Monocyte_121015_Long_C1D8","Monocyte_121015_Long_Post","Monocyte_118877_Long_Pre","Monocyte_118877_Long_C1D8","Monocyte_118877_Long_Post")]

library("circlize") ## For color options
DESeq2_NormCounts_broad_PBMC_GOI2 <- DESeq2_NormCounts_broad_PBMC_GOI2[order(DESeq2_NormCounts_broad_PBMC_GOI2$Gene),]
name <- DESeq2_NormCounts_broad_PBMC_GOI2[,c("Gene")]
df.OG2 <- data.matrix(log2(DESeq2_NormCounts_broad_PBMC_GOI2[,c(2:18)]+0.01))
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(DESeq2_NormCounts_broad_PBMC_GOI2[,c(2:18)])
pheatmap(df.OG2,color = colorRampPalette(c("#f0f0f0",rev(brewer.pal(n = 11, name = "Spectral"))))(100),breaks = c(seq(-0.1,15,by=0.151)),show_rownames = T,cluster_rows =F, cluster_cols=F)

## Generate Z-score for TPM from all samples (USE THIS)###
library(matrixStats)
DESeq2_NormCounts_broad_PBMC_GOI21<-log2(DESeq2_NormCounts_broad_PBMC_GOI2[,c(2:18)]+0.01)
DESeq2_NormCounts_broad_PBMC_GOI21_Zscore<- (DESeq2_NormCounts_broad_PBMC_GOI21-rowMeans(DESeq2_NormCounts_broad_PBMC_GOI21))/(rowSds(as.matrix(DESeq2_NormCounts_broad_PBMC_GOI21)))[row(DESeq2_NormCounts_broad_PBMC_GOI21)]
name <- DESeq2_NormCounts_broad_PBMC_GOI2[,c("Gene")]
df.OG2 <- data.matrix(DESeq2_NormCounts_broad_PBMC_GOI21_Zscore)
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(DESeq2_NormCounts_broad_PBMC_GOI2[,c(2:18)])
pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = T,cluster_rows =F, cluster_cols=F)






######Viral Mimicry  pathway broad cell annotaiton#####
GOIgenes <- read.table("~/jones-secondary/projects/JJANG/projects/CITEseq/Bulk_CD34_totalRNA/GSEA gene set lists/Minmin_Viral_mimicry_genelist.txt", sep = "\t", header = F, stringsAsFactors = F)

DESeq2_NormCounts_broad_PBMC_GOI <- DESeq2_NormCounts_broad_PBMC[DESeq2_NormCounts_broad_PBMC$Gene %in% GOIgenes$V1,]

DESeq2_NormCounts_broad_PBMC_GOI1 <- DESeq2_NormCounts_broad_PBMC_GOI
DESeq2_NormCounts_broad_PBMC_GOI1[is.na(DESeq2_NormCounts_broad_PBMC_GOI1)] <- -0.1
DESeq2_NormCounts_broad_PBMC_GOI1 <- DESeq2_NormCounts_broad_PBMC_GOI1[,c("Gene","CD34_24780_Short_Pre","CD34_24780_Short_Post","CD34_24780_Short_C2D8", "CD34_121474_Long_Pre","CD34_121474_Long_C1D8","CD34_121474_Long_Post","CD34_121015_Long_Pre","CD34_121015_Long_Post","CD34_118877_Long_Post","Monocyte_121474_Long_Pre","Monocyte_121474_Long_C1D8","Monocyte_121015_Long_Pre","Monocyte_121015_Long_C1D8","Monocyte_121015_Long_Post","Monocyte_118877_Long_Pre","Monocyte_118877_Long_C1D8","Monocyte_118877_Long_Post","B.cell_24780_Short_Pre","B.cell_24780_Short_Post","B.cell_24780_Short_C2D8", "B.cell_121474_Long_Pre","B.cell_121474_Long_C1D8","B.cell_121474_Long_Post", "B.cell_121015_Long_Pre","B.cell_121015_Long_C1D8","B.cell_121015_Long_Post", "B.cell_118877_Long_Pre","B.cell_118877_Long_C1D8","B.cell_118877_Long_Post","CD4_24780_Short_Pre","CD4_24780_Short_Post","CD4_24780_Short_C2D8", "CD4_121474_Long_Pre","CD4_121474_Long_C1D8","CD4_121474_Long_Post", "CD4_121015_Long_Pre","CD4_121015_Long_C1D8","CD4_121015_Long_Post", "CD4_118877_Long_Pre","CD4_118877_Long_C1D8","CD4_118877_Long_Post","CD8_24780_Short_Pre","CD8_24780_Short_Post","CD8_24780_Short_C2D8", "CD8_121474_Long_Pre","CD8_121474_Long_C1D8","CD8_121474_Long_Post", "CD8_121015_Long_Pre","CD8_121015_Long_C1D8","CD8_121015_Long_Post", "CD8_118877_Long_Pre","CD8_118877_Long_C1D8","CD8_118877_Long_Post","NK_24780_Short_Pre","NK_24780_Short_Post","NK_24780_Short_C2D8", "NK_121474_Long_Pre","NK_121474_Long_C1D8","NK_121474_Long_Post", "NK_121015_Long_Pre","NK_121015_Long_C1D8","NK_121015_Long_Post", "NK_118877_Long_Pre","NK_118877_Long_C1D8","NK_118877_Long_Post","gdT_24780_Short_Pre","gdT_24780_Short_Post","gdT_24780_Short_C2D8", "gdT_121474_Long_Pre","gdT_121474_Long_C1D8", "gdT_121015_Long_Pre","gdT_121015_Long_C1D8","gdT_121015_Long_Post", "gdT_118877_Long_Pre","gdT_118877_Long_C1D8")]
write.table(DESeq2_NormCounts_broad_PBMC_GOI1[,2:76]+0.01, "ViralMimicry_PseudoBUlk_PBMC_normCounts.txt", sep = "\t", col.names = F, row.names = F, quote = F)
#awk -v OFS="\t" '{print $1/$1,$2/$1,$3/$2,$4/$4,$5/$4,$6/$4,$7/$7,$8/$7,$9/$9,$10/$10,$11/$10,$12/$12,$13/$12,$14/$12,$15/$15,$16/$15,$17/$15,$18/$18,$19/$18,$20/$19,$21/$21,$22/$21,$23/$21,$24/$24,$25/$24,$26/$24,$27/$27,$28/$27,$29/$27,$30/$30,$31/$30,$32/$31,$33/$33,$34/$33,$35/$33, $36/$36,$37/$36,$38/$36, $39/$39,$40/$39,$41/$39, $42/$42,$43/$42,$44/$43,$45/$45,$46/$45,$47/$45, $48/$48, $49/$48, $50/$48, $51/$51, $52/$51, $53/$51, $54/$54,$55/$54,$56/$55, $57/$57, $58/$57,$59/$57, $60/$60, $61/$60,$62/$60, $63/$63, $64/$63,$65/$63, $66/$66, $67/$66, $68/$67,$69/$69, $70/$69, $71/$71, $72/$71, $73/$71,$74/$74, $75/$74 }' ViralMimicry_PseudoBUlk_PBMC_normCounts.txt > ViralMimicry_PseudoBUlk_PBMC_FC_analysis.txt

DESeq2_NormCounts_broad_PBMC_GOI2 <- DESeq2_NormCounts_broad_PBMC_GOI1[,c("Gene","CD34_24780_Short_Pre","CD34_24780_Short_Post","CD34_24780_Short_C2D8", "CD34_121474_Long_Pre","CD34_121474_Long_C1D8","CD34_121474_Long_Post","CD34_121015_Long_Pre","CD34_121015_Long_Post","CD34_118877_Long_Post")]

library("circlize") ## For color options
name <- DESeq2_NormCounts_broad_PBMC_GOI2[,c("Gene")]
df.OG2 <- data.matrix(log2(DESeq2_NormCounts_broad_PBMC_GOI2[,c(2:10)]+0.01))
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(DESeq2_NormCounts_broad_PBMC_GOI2[,c(2:10)])
pheatmap(df.OG2,color = colorRampPalette(c("#f0f0f0",rev(brewer.pal(n = 11, name = "Spectral"))))(100),breaks = c(seq(-0.1,15,by=0.151)),show_rownames = T,cluster_rows =F, cluster_cols=F)

## Generate Z-score for TPM from all samples (USE THIS)###
library(matrixStats)
DESeq2_NormCounts_broad_PBMC_GOI21<-log2(DESeq2_NormCounts_broad_PBMC_GOI2[,c(2:10)]+0.01)
DESeq2_NormCounts_broad_PBMC_GOI21_Zscore<- (DESeq2_NormCounts_broad_PBMC_GOI21-rowMeans(DESeq2_NormCounts_broad_PBMC_GOI21))/(rowSds(as.matrix(DESeq2_NormCounts_broad_PBMC_GOI21)))[row(DESeq2_NormCounts_broad_PBMC_GOI21)]
name <- DESeq2_NormCounts_broad_PBMC_GOI2[,c("Gene")]
df.OG2 <- data.matrix(DESeq2_NormCounts_broad_PBMC_GOI21_Zscore)
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(DESeq2_NormCounts_broad_PBMC_GOI2[,c(2:10)])
pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = T,cluster_rows =F, cluster_cols=F)


#### fold change heatmap
VM_FC <- read.table("ViralMimicry_PseudoBUlk_PBMC_FC_analysis.txt", header = F, sep ="\t", stringsAsFactors = F)
colnames(VM_FC) <- colnames(DESeq2_NormCounts_broad_PBMC_GOI1[,2:76])
VM_FC <- cbind(DESeq2_NormCounts_broad_PBMC_GOI1[,c("Gene")],VM_FC)

name <- VM_FC[,1]
df.OG2 <- data.matrix(log2(VM_FC[,2:76]))
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(VM_FC[,2:76])
pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PRGn")))(100),breaks = c(seq(-2,2,by=0.04)),clustering_method="ward.D2",show_rownames = T,cluster_rows =F, cluster_cols=F)







##broad anno (BOXPLOTS)
GOIgenes <- c("S100A8","S100A9") #Monocyte


DESeq2_NormCounts_broad_PBMC_GOI <- DESeq2_NormCounts_broad_PBMC[DESeq2_NormCounts_broad_PBMC$Gene %in% GOIgenes,]

m_broad_PBMC_GOI <- reshape2::melt(DESeq2_NormCounts_broad_PBMC_GOI, id = "Gene")
m_broad_PBMC_GOI[is.na(m_broad_PBMC_GOI)] <- 0
m_broad_PBMC_GOI$Anno <- apply(m_broad_PBMC_GOI, 1, function(x) {unlist(strsplit(x[2], "_"))[1]})
m_broad_PBMC_GOI$Patient <- apply(m_broad_PBMC_GOI, 1, function(x) {unlist(strsplit(x[2], "_"))[2]})
m_broad_PBMC_GOI$Survival <- apply(m_broad_PBMC_GOI, 1, function(x) {unlist(strsplit(x[2], "_"))[3]})
m_broad_PBMC_GOI$Time <- apply(m_broad_PBMC_GOI, 1, function(x) {unlist(strsplit(x[2], "_"))[4]})
m_broad_PBMC_GOI$SurvivalTime <- paste(m_broad_PBMC_GOI$Survival,m_broad_PBMC_GOI$Time, sep = "-")
m_broad_PBMC_GOI$PatientSurvivalTime <- paste(m_broad_PBMC_GOI$Patient,m_broad_PBMC_GOI$Survival,m_broad_PBMC_GOI$Time, sep = "-")


m_broad_PBMC_GOI$Anno <- factor(m_broad_PBMC_GOI$Anno, levels = c("CD34","Monocyte","B.cell","CD4","CD8","NK","gdT"))
m_broad_PBMC_GOI$SurvivalTime <- factor(m_broad_PBMC_GOI$SurvivalTime, levels = c("Short-Pre","Short-Post","Short-C2D8","Long-Pre","Long-C1D8","Long-Post"))

m_broad_PBMC_GOI$Gene <- factor(m_broad_PBMC_GOI$Gene, levels = GOIgenes)

p<-ggplot(m_broad_PBMC_GOI[m_broad_PBMC_GOI$Anno %in% c("CD34","Monocyte"),])+geom_boxplot(aes(Anno,value, color=SurvivalTime),outlier.shape = NA,width = 0.8)+geom_point(aes(Anno,value, group = SurvivalTime, color=SurvivalTime),size=1.5,position = position_dodge(width=0.8))+ggtitle("Immune Frequency Dynamics (Broad)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Anno", y = "normCounts")+scale_fill_manual(values =c("#d4d0d1","#A1A1A1","#5f5e60","#eecbff","#A24886","#8247a3"))+scale_color_manual(values =c("#d4d0d1","#A1A1A1","#5f5e60","#eecbff","#A24886","#8247a3"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
p+ facet_wrap( ~ Gene, scales="free",ncol=3)


##specific anno (BOXPLOTS)
DESeq2_NormCounts_specific_PBMC_GOI <- DESeq2_NormCounts_specific_PBMC[DESeq2_NormCounts_specific_PBMC$Gene %in% GOIgenes,]

m_specific_PBMC_GOI <- reshape2::melt(DESeq2_NormCounts_specific_PBMC_GOI, id = "Gene")
m_specific_PBMC_GOI[is.na(m_specific_PBMC_GOI)] <- 0
m_specific_PBMC_GOI$Anno <- apply(m_specific_PBMC_GOI, 1, function(x) {unlist(strsplit(x[2], "_"))[1]})
m_specific_PBMC_GOI$Patient <- apply(m_specific_PBMC_GOI, 1, function(x) {unlist(strsplit(x[2], "_"))[2]})
m_specific_PBMC_GOI$Survival <- apply(m_specific_PBMC_GOI, 1, function(x) {unlist(strsplit(x[2], "_"))[3]})
m_specific_PBMC_GOI$Time <- apply(m_specific_PBMC_GOI, 1, function(x) {unlist(strsplit(x[2], "_"))[4]})
m_specific_PBMC_GOI$SurvivalTime <- paste(m_specific_PBMC_GOI$Survival,m_specific_PBMC_GOI$Time, sep = "-")
m_specific_PBMC_GOI$PatientSurvivalTime <- paste(m_specific_PBMC_GOI$Patient,m_specific_PBMC_GOI$Survival,m_specific_PBMC_GOI$Time, sep = "-")


m_specific_PBMC_GOI$Anno <- factor(m_specific_PBMC_GOI$Anno, levels = c("EPC","MEP","HSC","LMPP","GMP","Classical.Mono","Intermediate.Mono","Non-classical.Mono","mDC","pDC","Naive.B","Memory.B","Plasmablast","CD56dim.NK","CD56bright.NK","Vd1.gdT","Vd2.gdT","CD4.Naive","CD4.TCM","CD4.TEM","CD4.CTL","CD4.Treg","CD8.Naive","CD8.TCM","CD8.TEM","CD8.TEMRA","MAIT",'Stromal',"Megakaryocyte","Neutrophil","Differentiating.Stem.Cell"))
m_specific_PBMC_GOI$SurvivalTime <- factor(m_specific_PBMC_GOI$SurvivalTime, levels = c("Short-Pre","Short-Post","Short-C2D8","Long-Pre","Long-C1D8","Long-Post"))

GOIgenes <- c("IFNG")

m_specific_PBMC_GOI$Gene <- factor(m_specific_PBMC_GOI$Gene, levels = GOIgenes)


p<-ggplot(m_specific_PBMC_GOI)+geom_boxplot(aes(Anno,value, color=SurvivalTime),outlier.shape = NA,width = 0.8)+geom_point(aes(Anno,value, group = SurvivalTime, color=SurvivalTime),size=1.5,position = position_dodge(width=0.8))+ggtitle("Immune Frequency Dynamics (Broad)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Anno", y = "normCounts")+scale_fill_manual(values =c("#d4d0d1","#A1A1A1","#5f5e60","#eecbff","#A24886","#8247a3"))+scale_color_manual(values =c("#d4d0d1","#A1A1A1","#5f5e60","#eecbff","#A24886","#8247a3"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
p+ facet_wrap( ~ Gene, scales="free",ncol=3)





#### PLOT EXHAUSTION MARKERS (ADT) FOR PBMC cells ####

DefaultAssay(BM.combined.sct.PBMC) <-"ADT2"
Idents(BM.combined.sct.PBMC) <- "manualAnno_specific2"
BM.combined.sct.PBMC$manualAnno_specific2 <- factor(BM.combined.sct.PBMC$manualAnno_specific2, levels = c("EPC","MEP","HSC","LMPP","GMP","Classical Mono","Intermediate Mono","Non-classical Mono","mDC","pDC","Naive B","Memory B","Plasmablast","CD56dim NK","CD56bright NK","Vd1 gdT","Vd2 gdT","CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL","CD4 Treg","CD8 Naive","CD8 TCM","CD8 TEM","CD8 TEMRA","MAIT",'Stromal', "Megakaryocyte","Neutrophil","Differentiating Stem Cell"))
BM.combined.sct.PBMC$TimeSurvival <- factor(BM.combined.sct.PBMC$TimeSurvival, levels = c("Short_Pre","Short_Post","Short_C2D8","Long_Pre","Long_C1D8","Long_Post"))
levels(BM.combined.sct.PBMC) <- c("EPC","MEP","HSC","LMPP","GMP","Classical Mono","Intermediate Mono","Non-classical Mono","mDC","pDC","Naive B","Memory B","Plasmablast","CD56dim NK","CD56bright NK","Vd1 gdT","Vd2 gdT","CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL","CD4 Treg","CD8 Naive","CD8 TCM","CD8 TEM","CD8 TEMRA","MAIT",'Stromal', "Megakaryocyte","Neutrophil","Differentiating Stem Cell")
  
aa<-VlnPlot(object = BM.combined.sct.PBMC, layer = "data", features = c("CD279"),idents = c("HSC","LMPP","Naive B","Memory B","CD56dim NK","Vd1 gdT","Vd2 gdT","CD4 Naive","CD4 TCM","CD8 Naive","CD8 TEM","CD8 TEMRA"), cols = c("#d4d0d1","#A1A1A1","#5f5e60","#eecbff","#A24886","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
bb<-VlnPlot(object = BM.combined.sct.PBMC, layer = "data", features = c("TIGIT-VSTM3"),idents = c("HSC","LMPP","Naive B","Memory B","CD56dim NK","Vd1 gdT","Vd2 gdT","CD4 Naive","CD4 TCM","CD8 Naive","CD8 TEM","CD8 TEMRA"), cols = c("#d4d0d1","#A1A1A1","#5f5e60","#eecbff","#A24886","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
cc<-VlnPlot(object = BM.combined.sct.PBMC, layer = "data", features = c("CD223-LAG3"),idents = c("HSC","LMPP","Naive B","Memory B","CD56dim NK","Vd1 gdT","Vd2 gdT","CD4 Naive","CD4 TCM","CD8 Naive","CD8 TEM","CD8 TEMRA"), cols = c("#d4d0d1","#A1A1A1","#5f5e60","#eecbff","#A24886","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
dd<-VlnPlot(object = BM.combined.sct.PBMC, layer = "data", features = c("CD366-Tim3"),idents = c("HSC","LMPP","Naive B","Memory B","CD56dim NK","Vd1 gdT","Vd2 gdT","CD4 Naive","CD4 TCM","CD8 Naive","CD8 TEM","CD8 TEMRA"), cols = c("#d4d0d1","#A1A1A1","#5f5e60","#eecbff","#A24886","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')

multiplot(aa,bb,cc,dd,cols=1)


#Perform bone marrow specific with patients profiled also with PBMC

DefaultAssay(BM.combined.sct.BM) <-"ADT2"
Idents(BM.combined.sct.BM) <- "manualAnno_specific2"
BM.combined.sct.BM$manualAnno_specific2 <- factor(BM.combined.sct.BM$manualAnno_specific2, levels = c("EPC","MEP","HSC","LMPP","GMP","Classical Mono","Intermediate Mono","Non-classical Mono","mDC","pDC","Naive B","Memory B","Plasmablast","CD56dim NK","CD56bright NK","Vd1 gdT","Vd2 gdT","CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL","CD4 Treg","CD8 Naive","CD8 TCM","CD8 TEM","CD8 TEMRA","MAIT",'Stromal', "Megakaryocyte","Neutrophil","Differentiating Stem Cell"))
BM.combined.sct.BM$TimeSurvival <- factor(BM.combined.sct.BM$TimeSurvival, levels = c("Short_Pre","Short_Post","Short_C2D8","Long_Pre","Long_C1D8","Long_Post"))
levels(BM.combined.sct.BM) <- c("EPC","MEP","HSC","LMPP","GMP","Classical Mono","Intermediate Mono","Non-classical Mono","mDC","pDC","Naive B","Memory B","Plasmablast","CD56dim NK","CD56bright NK","Vd1 gdT","Vd2 gdT","CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL","CD4 Treg","CD8 Naive","CD8 TCM","CD8 TEM","CD8 TEMRA","MAIT",'Stromal', "Megakaryocyte","Neutrophil","Differentiating Stem Cell")

Idents(BM.combined.sct.BM) <- "Patient"
BM.combined.sct.BM_subset <- subset(x = BM.combined.sct.BM, idents = c("24780","121474","121015","118877"))


Idents(BM.combined.sct.BM_subset) <- "manualAnno_specific2"
aa<-VlnPlot(object = BM.combined.sct.BM_subset, layer = "data", features = c("CD279"),idents = c("HSC","LMPP","Naive B","Memory B","CD56dim NK","Vd1 gdT","Vd2 gdT","CD4 Naive","CD4 TCM","CD8 Naive","CD8 TEM","CD8 TEMRA"), cols = c("#d4d0d1","#A1A1A1","#5f5e60","#eecbff","#A24886","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
bb<-VlnPlot(object = BM.combined.sct.BM_subset, layer = "data", features = c("TIGIT-VSTM3"),idents = c("HSC","LMPP","Naive B","Memory B","CD56dim NK","Vd1 gdT","Vd2 gdT","CD4 Naive","CD4 TCM","CD8 Naive","CD8 TEM","CD8 TEMRA"), cols = c("#d4d0d1","#A1A1A1","#5f5e60","#eecbff","#A24886","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
cc<-VlnPlot(object = BM.combined.sct.BM_subset, layer = "data", features = c("CD223-LAG3"),idents = c("HSC","LMPP","Naive B","Memory B","CD56dim NK","Vd1 gdT","Vd2 gdT","CD4 Naive","CD4 TCM","CD8 Naive","CD8 TEM","CD8 TEMRA"), cols = c("#d4d0d1","#A1A1A1","#5f5e60","#eecbff","#A24886","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')
dd<-VlnPlot(object = BM.combined.sct.BM_subset, layer = "data", features = c("CD366-Tim3"),idents = c("HSC","LMPP","Naive B","Memory B","CD56dim NK","Vd1 gdT","Vd2 gdT","CD4 Naive","CD4 TCM","CD8 Naive","CD8 TEM","CD8 TEMRA"), cols = c("#d4d0d1","#A1A1A1","#5f5e60","#eecbff","#A24886","#8247a3"),split.by = 'TimeSurvival', ncol=1, pt.size=0,alpha=0.1)+ geom_boxplot(position=position_dodge(0.9),outlier.shape = NA,alpha =0.1) + theme(legend.position = 'none')

multiplot(aa,bb,cc,dd,cols=1)
