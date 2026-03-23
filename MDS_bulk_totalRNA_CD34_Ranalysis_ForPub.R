library(tibble)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library("circlize") ## For color options
library(reshape2)
library(matrixStats)
library(pheatmap)
library(dplyr)
library(DESeq2)
library("viridis")
library(EnhancedVolcano)
library(cowplot)
library(fgsea)
library(tidyr)


mypalette <- brewer.pal(12,"Paired")
mypalette2 <- brewer.pal(12,"Set3")
mypalette3 <- brewer.pal(8,"Pastel2")
mypalette4 <- brewer.pal(8,"Set2")
mypalette5 <- brewer.pal(8,"Pastel1")
mypalette6 <- brewer.pal(9,"Greys")[2:9]

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





####### QC DATA (Long = >15 MOS, Short = <15 months)    ###############
setwd("~/Bulk_CD34_totalRNA_forPub")
meta <- read.table("MDS_meta_final_forPub.txt", sep = "\t", header = T, stringsAsFactors = F)
meta$Time2 <- factor(meta$Time2, levels = c("Scr","PC2","PC4"))
meta$Survival<- factor(meta$Survival, levels = c("Healthy","Short","Long"))
meta$SurvivalTime<- factor(meta$SurvivalTime, levels = c("HealthySc","ShortSc","Short2","Short4","LongSc","Long2","Long4"))
meta<-meta[order(meta$Survival,meta$MOS,meta$SubjectID,meta$Time),]
meta$Sample<- factor(meta$Sample, levels = c(as.character(meta$Sample)))

meta_paired <- meta[meta$SubjectID %in% names(which(table(as.character(meta$SubjectID)) > 1)),]
meta_paired$SurvivalTime<- factor(meta_paired$SurvivalTime, levels = c("HealthySc","ShortSc","Short2","Short4","LongSc","Long2","Long4"))




############## DESeq 2 Analysis for PERFORM edgeR TMM NORmalization ON RAW COUNTS and for DEG calling (USE THIS)###########
Ensb_genelist<-read.table("ENSEMBL_gene_transcript_ID_Genename.txt", sep = "\t", header =T, stringsAsFactors = F)
colnames(Ensb_genelist)[4]<- "GeneID"

countData <- as.matrix(read.csv("gene_count_matrix_forpub.csv", row.names="gene_id"))
colData <- meta
rownames(colData) <-as.character(colData$Sample)
countData <- countData[,colnames(countData) %in% as.character(colData$Sample)]
all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))
colData$SurvivalTime <- as.character(colData$SurvivalTime)
colData$Batch <- as.factor(colData$Batch)



dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~SurvivalTime)
keep <- rowSums(counts(dds)) >= 10 #filter lowly expressed genes
dds <- dds[keep,]
dds <- DESeq(dds)



cdds<-data.frame(counts(dds, normalized = TRUE))
cdds$GeneID <- rownames(cdds)
cdds$Gene.name <- apply(cdds, 1, function(x) {unlist(strsplit(x[51],"\\|"))[2]})
cdds<-cdds[,c(51,52,1:50)]
cdds$GeneID <- apply(cdds, 1, function(x) {unlist(strsplit(x[1],"\\|"))[1]})


### Need to analyze samples with paired timepoints only
countData_Short2 <- as.matrix(read.csv("gene_count_matrix_forpub.csv", row.names="gene_id"))
colData_Short2 <- meta_paired[meta_paired$Survival == "Short" & (meta_paired$Time2 == "Scr" | meta_paired$Time2 == "PC2"),]
colData_Short2 <- colData_Short2[colData_Short2$SubjectID %in% c("121589","113446","116467"),]
rownames(colData_Short2 ) <-as.character(colData_Short2$Sample)
countData_Short2  <- countData_Short2[,colnames(countData_Short2) %in% as.character(colData_Short2$Sample)]
all(rownames(colData_Short2) %in% colnames(countData_Short2))
countData_Short2 <- countData_Short2[, rownames(colData_Short2)]
all(rownames(colData_Short2) == colnames(countData_Short2))
colData_Short2$SurvivalTime <- as.character(colData_Short2$SurvivalTime)
colData_Short2$Batch <- as.character(colData_Short2$Batch)

countData_Short4 <- as.matrix(read.csv("gene_count_matrix_forpub.csv", row.names="gene_id"))
colData_Short4 <-  meta_paired[meta_paired$Survival == "Short" & (meta_paired$Time2 == "Scr" | meta_paired$Time2 == "PC4"),]
colData_Short4 <- colData_Short4[colData_Short4$SubjectID %in% c("116920","113446","116467"),]
rownames(colData_Short4 ) <-as.character(colData_Short4$Sample)
countData_Short4  <- countData_Short4[,colnames(countData_Short4) %in% as.character(colData_Short4$Sample)]
all(rownames(colData_Short4) %in% colnames(countData_Short4))
countData_Short4 <- countData_Short4[, rownames(colData_Short4)]
all(rownames(colData_Short4) == colnames(countData_Short4))
colData_Short4$SurvivalTime <- as.character(colData_Short4$SurvivalTime)
colData_Short4$Batch <- as.character(colData_Short4$Batch)

countData_Long2 <- as.matrix(read.csv("gene_count_matrix_forpub.csv", row.names="gene_id"))
colData_Long2 <- meta_paired[meta_paired$Survival == "Long" & (meta_paired$Time2 == "Scr" | meta_paired$Time2 == "PC2"),]
colData_Long2 <- colData_Long2[colData_Long2$SubjectID %in% c("121474","122739","117494","121535","119837","121015","114152","119325","118877","115585"),]
rownames(colData_Long2 ) <-as.character(colData_Long2$Sample)
countData_Long2  <- countData_Long2[,colnames(countData_Long2) %in% as.character(colData_Long2$Sample)]
all(rownames(colData_Long2) %in% colnames(countData_Long2))
countData_Long2 <- countData_Long2[, rownames(colData_Long2)]
all(rownames(colData_Long2) == colnames(countData_Long2))
colData_Long2$SurvivalTime <- as.character(colData_Long2$SurvivalTime)
colData_Long2$Batch <- as.character(colData_Long2$Batch)

countData_Long4 <- as.matrix(read.csv("gene_count_matrix_forpub.csv", row.names="gene_id"))
colData_Long4 <- meta_paired[meta_paired$Survival == "Long" & (meta_paired$Time2 == "Scr" | meta_paired$Time2 == "PC4"),]
colData_Long4 <- colData_Long4[colData_Long4$SubjectID %in% c("117494","121535","119837","121015","114152"),]
rownames(colData_Long4 ) <-as.character(colData_Long4$Sample)
countData_Long4  <- countData_Long4[,colnames(countData_Long4) %in% as.character(colData_Long4$Sample)]
all(rownames(colData_Long4) %in% colnames(countData_Long4))
countData_Long4 <- countData_Long4[, rownames(colData_Long4)]
all(rownames(colData_Long4) == colnames(countData_Long4))
colData_Long4$SurvivalTime <- as.character(colData_Long4$SurvivalTime)
colData_Long4$Batch <- as.character(colData_Long4$Batch)


dds_Short2 <- DESeqDataSetFromMatrix(countData = countData_Short2, colData = colData_Short2, design = ~ SurvivalTime)
keep <- rowSums(counts(dds_Short2)) >= 10 #filter lowly expressed genes
dds_Short2 <- dds_Short2[keep,]
dds_Short2 <- DESeq(dds_Short2)

dds_Short4 <- DESeqDataSetFromMatrix(countData = countData_Short4, colData = colData_Short4, design = ~SurvivalTime)
keep <- rowSums(counts(dds_Short4)) >= 10 #filter lowly expressed genes
dds_Short4 <- dds_Short4[keep,]
dds_Short4 <- DESeq(dds_Short4)

dds_Long2 <- DESeqDataSetFromMatrix(countData = countData_Long2, colData = colData_Long2, design = ~SurvivalTime)
keep <- rowSums(counts(dds_Long2)) >= 10 #filter lowly expressed genes
dds_Long2 <- dds_Long2[keep,]
dds_Long2 <- DESeq(dds_Long2)

dds_Long4 <- DESeqDataSetFromMatrix(countData = countData_Long4, colData = colData_Long4, design = ~SurvivalTime)
keep <- rowSums(counts(dds_Long4)) >= 10 #filter lowly expressed genes
dds_Long4 <- dds_Long4[keep,]
dds_Long4 <- DESeq(dds_Long4)


vsd <- vst(dds, blind=FALSE)
mat <- assay(vsd)
mat <- limma::removeBatchEffect(mat, vsd$Batch)
assay(vsd) <- mat
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( c(brewer.pal(5, "Spectral")) )(255)
d<-pheatmap(sampleDistMatrix,
            clustering_distance_rows=sampleDists,
            clustering_distance_cols=sampleDists,
            col=colors)
meta_cluster <- meta[meta$Sample %in% d$tree_col$labels,]
meta_cluster$SubjectID<-factor(meta_cluster$SubjectID, levels =c(unique(as.character(meta$SubjectID))))
meta_cluster<-meta_cluster[d$tree_col$order,]
meta_cluster$Sample<- factor(meta_cluster$Sample, levels = c(as.character(meta_cluster$Sample)))
meta_cluster$SurvivalTime<- factor(meta_cluster$SurvivalTime, levels = c("HealthySc","ShortSc", "Short2","Short4","LongSc","Long2","Long4"))

d3<-ggplot(meta_cluster, aes(y = 1, x = Sample, fill = as.character(SubjectID)))+geom_tile()+scale_fill_manual(values = c(mypalette,mypalette3,"pink",rep("#A74E60",8)))+scale_color_manual(values=c(mypalette,mypalette3,"black"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
d7<-ggplot(meta_cluster, aes(y = 1, x = Sample, fill = SurvivalTime))+scale_fill_manual(values=c("#A74E60", "#c0c5ce","#65737e","black" ,"#deacde","#be75be","#872e87"))+geom_tile()+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))

multiplot(d3,d7,cols=1)


a<-plotPCA(vsd, intgroup=c("SubjectID"))+scale_color_manual(values=c(mypalette,mypalette3,"pink",rep("#A74E60",8)))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
b<-plotPCA(vsd, intgroup=c("SurvivalTime"))+scale_color_manual(values=c("#A74E60","#65737e","black","#c0c5ce" ,"#be75be","#872e87","#deacde"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))

multiplot(a,b,cols=1)




###################### DIFFERENTIAL GENE CALLING (DESeq2) ###########################
#perform DESeq on all samples
res_Short_scr<- results(dds,name="SurvivalTime_ShortSc_vs_LongSc", contrast=c("SurvivalTime","ShortSc","LongSc"))
res_Short_scr<- res_Short_scr[order(-res_Short_scr$baseMean),]

DEG_res_Short_scr<- data.frame(res_Short_scr)
DEG_res_Short_scr$geneinfo <- rownames(DEG_res_Short_scr)
DEG_res_Short_scr$gene <- apply(DEG_res_Short_scr, 1, function(x) {unlist(strsplit(x[7],"\\|"))[2]})

DEG_res_Short_scr_0.1<- DEG_res_Short_scr[which(DEG_res_Short_scr$padj < 0.1), ]
DEG_res_Short_scr_0.1$Comparison <- "Short Survivor Screening vs. Long Survivor Screening"
DEG_res_Short_scr_0.1 <- DEG_res_Short_scr_0.1[order(DEG_res_Short_scr_0.1$log2FoldChange),]


v1 <- EnhancedVolcano(DEG_res_Short_scr,
                lab = DEG_res_Short_scr$gene,
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'P.adj'),
                pCutoff = 0.1,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 5.0,
                colAlpha = 0.8,
                labFace = 'bold',
                 legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                col=c('grey', 'grey', 'grey', 'red3'),
                colConnectors = 'black',
                drawConnectors = F)




res_Short_scr_2_paired<- results(dds_Short2 ,name="SurvivalTime_ShortSc_vs_Short2", contrast=c("SurvivalTime","ShortSc","Short2"))
res_Short_scr_2_paired<- res_Short_scr_2_paired[order(-res_Short_scr_2_paired$baseMean),]

DEG_res_Short_scr_2_paired<- data.frame(res_Short_scr_2_paired)
DEG_res_Short_scr_2_paired$geneinfo <- rownames(DEG_res_Short_scr_2_paired)
DEG_res_Short_scr_2_paired$gene <- apply(DEG_res_Short_scr_2_paired, 1, function(x) {unlist(strsplit(x[7],"\\|"))[2]})
DEG_res_Short_scr_2_paired[is.na(DEG_res_Short_scr_2_paired)] <- 1

DEG_res_Short_scr_2_paired_0.1<- DEG_res_Short_scr_2_paired[which(DEG_res_Short_scr_2_paired$padj < 0.1), ]
DEG_res_Short_scr_2_paired_0.1$Comparison <- "Short Survivor Screening vs. Short Survivor PC2"
DEG_res_Short_scr_2_paired_0.1 <- DEG_res_Short_scr_2_paired_0.1[order(DEG_res_Short_scr_2_paired_0.1$log2FoldChange),]


v2 <- EnhancedVolcano(DEG_res_Short_scr_2_paired,
                lab = DEG_res_Short_scr_2_paired$gene,
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'P.adj'),
                pCutoff = 0.1,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 5.0,
                colAlpha = 0.8,
                labFace = 'bold',
                 legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                col=c('grey', 'grey', 'grey', 'red3'),
                colConnectors = 'black',
                drawConnectors = F)





res_Short_scr_4_paired<- results(dds_Short4,name="SurvivalTime_ShortSc_vs_Short4", contrast=c("SurvivalTime","ShortSc","Short4"))
res_Short_scr_4_paired<- res_Short_scr_4_paired[order(-res_Short_scr_4_paired$baseMean),]

DEG_res_Short_scr_4_paired<- data.frame(res_Short_scr_4_paired)
DEG_res_Short_scr_4_paired$geneinfo <- rownames(DEG_res_Short_scr_4_paired)
DEG_res_Short_scr_4_paired$gene <- apply(DEG_res_Short_scr_4_paired, 1, function(x) {unlist(strsplit(x[7],"\\|"))[2]})
DEG_res_Short_scr_4_paired[is.na(DEG_res_Short_scr_4_paired)] <- 1

DEG_res_Short_scr_4_paired_0.1<- DEG_res_Short_scr_4_paired[which(DEG_res_Short_scr_4_paired$padj < 0.1), ]
DEG_res_Short_scr_4_paired_0.1$Comparison <- "Short Survivor Screening vs. Short Survivor PC4"
DEG_res_Short_scr_4_paired_0.1 <- DEG_res_Short_scr_4_paired_0.1[order(DEG_res_Short_scr_4_paired_0.1$log2FoldChange),]



v3 <- EnhancedVolcano(DEG_res_Short_scr_4_paired,
                lab = DEG_res_Short_scr_4_paired$gene,
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'P.adj'),
                pCutoff = 0.1,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 5.0,
                colAlpha = 0.8,
                labFace = 'bold',
                 legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                col=c('grey', 'grey', 'grey', 'red3'),
                colConnectors = 'black',
                drawConnectors = F)



res_Long_scr_2_paired<- results(dds_Long2,name="SurvivalTime_LongSc_vs_Long2", contrast=c("SurvivalTime","LongSc","Long2"))
res_Long_scr_2_paired<- res_Long_scr_2_paired[order(-res_Long_scr_2_paired$baseMean),]

DEG_res_Long_scr_2_paired<- data.frame(res_Long_scr_2_paired)
DEG_res_Long_scr_2_paired$geneinfo <- rownames(DEG_res_Long_scr_2_paired)
DEG_res_Long_scr_2_paired$gene <- apply(DEG_res_Long_scr_2_paired, 1, function(x) {unlist(strsplit(x[7],"\\|"))[2]})
DEG_res_Long_scr_2_paired[is.na(DEG_res_Long_scr_2_paired)] <- 1

DEG_res_Long_scr_2_paired_0.1<- DEG_res_Long_scr_2_paired[which(DEG_res_Long_scr_2_paired$padj < 0.1), ]
DEG_res_Long_scr_2_paired_0.1$Comparison <- "Long Survivor Screening vs. Long Survivor PC2"
DEG_res_Long_scr_2_paired_0.1 <- DEG_res_Long_scr_2_paired_0.1[order(DEG_res_Long_scr_2_paired_0.1$log2FoldChange),]


v4 <- EnhancedVolcano(DEG_res_Long_scr_2_paired,
                lab = DEG_res_Long_scr_2_paired$gene,
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'P.adj'),
                pCutoff = 0.1,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 5.0,
                colAlpha = 0.8,
                labFace = 'bold',
                 legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                col=c('grey', 'grey', 'grey', 'red3'),
                colConnectors = 'black',
                drawConnectors = F)



res_Long_scr_4_paired<- results(dds_Long4,name="SurvivalTime_LongSc_vs_Long4", contrast=c("SurvivalTime","LongSc","Long4"))
res_Long_scr_4_paired<- res_Long_scr_4_paired[order(-res_Long_scr_4_paired$baseMean),]

DEG_res_Long_scr_4_paired<- data.frame(res_Long_scr_4_paired)
DEG_res_Long_scr_4_paired$geneinfo <- rownames(DEG_res_Long_scr_4_paired)
DEG_res_Long_scr_4_paired$gene <- apply(DEG_res_Long_scr_4_paired, 1, function(x) {unlist(strsplit(x[7],"\\|"))[2]})
DEG_res_Long_scr_4_paired[is.na(DEG_res_Long_scr_4_paired)] <- 1

DEG_res_Long_scr_4_paired_0.1<- DEG_res_Long_scr_4_paired[which(DEG_res_Long_scr_4_paired$padj < 0.1), ]
DEG_res_Long_scr_4_paired_0.1$Comparison <- "Long Survivor Screening vs. Long Survivor PC4"
DEG_res_Long_scr_4_paired_0.1 <- DEG_res_Long_scr_4_paired_0.1[order(DEG_res_Long_scr_4_paired_0.1$log2FoldChange),]

v5 <- EnhancedVolcano(DEG_res_Long_scr_4_paired,
                lab = DEG_res_Long_scr_4_paired$gene,
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'P.adj'),
                pCutoff = 0.1,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 5.0,
                colAlpha = 0.8,
                labFace = 'bold',
                 legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                col=c('grey', 'grey', 'grey', 'red3'),
                colConnectors = 'black',
                drawConnectors = F)




res_ShortvsLong_2<- results(dds,name="SurvivalTime_Short2_vs_Long2", contrast=c("SurvivalTime","Short2","Long2"))
res_ShortvsLong_2<- res_ShortvsLong_2[order(-res_ShortvsLong_2$baseMean),]

DEG_res_ShortvsLong_2<- data.frame(res_ShortvsLong_2)
DEG_res_ShortvsLong_2$geneinfo <- rownames(DEG_res_ShortvsLong_2)
DEG_res_ShortvsLong_2$gene <- apply(DEG_res_ShortvsLong_2, 1, function(x) {unlist(strsplit(x[7],"\\|"))[2]})
DEG_res_ShortvsLong_2[is.na(DEG_res_ShortvsLong_2)] <- 1

DEG_res_ShortvsLong_2_0.1<- DEG_res_ShortvsLong_2[which(DEG_res_ShortvsLong_2$padj < 0.1), ]
DEG_res_ShortvsLong_2_0.1$Comparison <- "Short Survivor PC2 vs. Long Survivor PC2"
DEG_res_ShortvsLong_2_0.1 <- DEG_res_ShortvsLong_2_0.1[order(DEG_res_ShortvsLong_2_0.1$log2FoldChange),]


v6 <- EnhancedVolcano(DEG_res_ShortvsLong_2,
                lab = DEG_res_ShortvsLong_2$gene,
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'P.adj'),
                pCutoff = 0.1,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 5.0,
                colAlpha = 0.8,
                labFace = 'bold',
                 legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                col=c('grey', 'grey', 'grey', 'red3'),
                colConnectors = 'black',
                drawConnectors = F)




res_ShortvsLong_4<- results(dds,name="Response4Time_Short4_vs_Long4", contrast=c("SurvivalTime","Short4","Long4"))
res_ShortvsLong_4<- res_ShortvsLong_4[order(-res_ShortvsLong_4$baseMean),]

DEG_res_ShortvsLong_4 <- data.frame(subset(res_ShortvsLong_4, padj < 0.1))
DEG_res_ShortvsLong_4<- data.frame(res_ShortvsLong_4)
DEG_res_ShortvsLong_4$geneinfo <- rownames(DEG_res_ShortvsLong_4)
DEG_res_ShortvsLong_4$gene <- apply(DEG_res_ShortvsLong_4, 1, function(x) {unlist(strsplit(x[7],"\\|"))[2]})
DEG_res_ShortvsLong_4[is.na(DEG_res_ShortvsLong_4)] <- 1

DEG_res_ShortvsLong_4_0.1<- DEG_res_ShortvsLong_4[which(DEG_res_ShortvsLong_4$padj < 0.1), ]
DEG_res_ShortvsLong_4_0.1$Comparison <- "Short Survivor PC4 vs. Long Survivor PC4"
DEG_res_ShortvsLong_4_0.1 <- DEG_res_ShortvsLong_4_0.1[order(DEG_res_ShortvsLong_4_0.1$log2FoldChange),]


v7 <- EnhancedVolcano(DEG_res_ShortvsLong_4,
                lab = DEG_res_ShortvsLong_4$gene,
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'P.adj'),
                pCutoff = 0.1,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 5.0,
                colAlpha = 0.8,
                labFace = 'bold',
                legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                col=c('grey', 'grey', 'grey', 'red3'),
                colConnectors = 'black',
                drawConnectors = F)

multiplot(v1,v1,v2,v3,v4,v5,v6,v7, cols = 2)


####Make bargraph with # of DEGs###
nrow(DEG_res_Short_scr_0.1[DEG_res_Short_scr_0.1$log2FoldChange > 1 ,]) #20
nrow(DEG_res_Short_scr_0.1[DEG_res_Short_scr_0.1$log2FoldChange < -1 ,]) #10
nrow(DEG_res_Short_scr_2_paired_0.1[DEG_res_Short_scr_2_paired_0.1$log2FoldChange > 1 ,]) #7
nrow(DEG_res_Short_scr_2_paired_0.1[DEG_res_Short_scr_2_paired_0.1$log2FoldChange < -1 ,]) #6
nrow(DEG_res_Short_scr_4_paired_0.1[DEG_res_Short_scr_4_paired_0.1$log2FoldChange > 1 ,]) #1
nrow(DEG_res_Short_scr_4_paired_0.1[DEG_res_Short_scr_4_paired_0.1$log2FoldChange < -1 ,]) #4
nrow(DEG_res_Long_scr_2_paired_0.1[DEG_res_Long_scr_2_paired_0.1$log2FoldChange > 1 ,]) #1
nrow(DEG_res_Long_scr_2_paired_0.1[DEG_res_Long_scr_2_paired_0.1$log2FoldChange < -1 ,]) #10
nrow(DEG_res_Long_scr_4_paired_0.1[DEG_res_Long_scr_4_paired_0.1$log2FoldChange > 1 ,]) #2
nrow(DEG_res_Long_scr_4_paired_0.1[DEG_res_Long_scr_4_paired_0.1$log2FoldChange < -1 ,]) #3
nrow(DEG_res_ShortvsLong_2_0.1[DEG_res_ShortvsLong_2_0.1$log2FoldChange > 1 ,]) #12
nrow(DEG_res_ShortvsLong_2_0.1[DEG_res_ShortvsLong_2_0.1$log2FoldChange < -1 ,]) #338
nrow(DEG_res_ShortvsLong_4_0.1[DEG_res_ShortvsLong_4_0.1$log2FoldChange > 1 ,]) #10
nrow(DEG_res_ShortvsLong_4_0.1[DEG_res_ShortvsLong_4_0.1$log2FoldChange < -1 ,]) #48

DEGbar <- data.frame(c("ShortScr_LongScr","ShortScr_LongScr","ShortScr_ShortPC2","ShortScr_ShortPC2","ShortScr_ShortPC4","ShortScr_ShortPC4","LongScr_LongPC2","LongScr_LongPC2","LongScr_LongPC4","LongScr_LongPC4","ShortPC2_LongPC2","ShortPC2_LongPC2","ShortPC4_LongPC4","ShortPC4_LongPC4"), c(rep(c("Up","Down"),7)),c(20,10,7,6,1,4,1,10,2,3,12,338,10,48))
colnames(DEGbar)<- c("Comparison","DEGdirection","GeneCount")
DEGbar$Comparison <- factor(DEGbar$Comparison, levels = rev(c("ShortScr_LongScr","ShortPC2_LongPC2","ShortPC4_LongPC4","ShortScr_ShortPC2","ShortScr_ShortPC4","LongScr_LongPC2","LongScr_LongPC4")))
DEGbar$DEGdirection <- factor(DEGbar$DEGdirection, levels = rev(c("Up","Down")))

dbar <- ggplot(DEGbar, aes(Comparison,GeneCount,fill=DEGdirection))+geom_bar(stat ="identity",position="dodge",colour = "black")+ggtitle("Number of DEGs")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  scale_y_continuous(limits = c(0,400))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+labs(x = "Comparisons", y = "# of DEGs")+scale_fill_manual(values = c("#667b68","#c99789"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
dbar + coord_flip()


#Make DEG supplementary table#
All_DEGs <- rbind(DEG_res_Short_scr_0.1,DEG_res_Short_scr_2_paired_0.1,DEG_res_Short_scr_4_paired_0.1,DEG_res_Long_scr_2_paired_0.1,DEG_res_Long_scr_4_paired_0.1,DEG_res_ShortvsLong_2_0.1,DEG_res_ShortvsLong_4_0.1)
All_DEGs <- All_DEGs[,c("gene","baseMean","log2FoldChange","padj","Comparison")]
All_DEGs_LFCcutoff1 <- All_DEGs[All_DEGs$log2FoldChange > 1 | All_DEGs$log2FoldChange < -1,] 

write.table(All_DEGs_LFCcutoff1, "All_DEGs_LFCcutoff1.txt", sep = "\t", row.names = F, col.names = T, quote = F )





#######################  GSEA Enrichment    ################################
Ensb_genelist<-read.table("../GSEA gene set lists/ENSEMBL_gene_transcript_ID_Genename.txt", sep = "\t", header =T, stringsAsFactors = F)
ens2symbol <- as_tibble(Ensb_genelist)
colnames(ens2symbol)[4] <- "GeneID"
pathways.hallmark <- gmtPathways("../GSEA gene set lists/h.all.v7.1.symbols.gmt")
pathways.hallmark <- gmtPathways("../GSEA gene set lists/c5.all.v7.1.symbols.gmt") # GO terms

fgsea_SENMAYO <- list(c("ACVR1B","ANG","ANGPT1","ANGPTL4","AREG","AXL","BEX3","BMP2","BMP6","C3","CCL1","CCL13","CCL16","CCL2","CCL20","CCL24","CCL26","CCL3","CCL3L1","CCL4","CCL5","CCL7","CCL8","CD55","CD9","CSF1","CSF2","CSF2RB","CST4","CTNNB1","CTSB","CXCL1","CXCL10","CXCL12","CXCL16","CXCL2","CXCL3","CXCL8","CXCR2","DKK1","EDN1","EGF","EGFR","EREG","ESM1","ETS2","FAS","FGF1","FGF2","FGF7","GDF15","GEM","GMFG","HGF","HMGB1","ICAM1","ICAM3","IGF1","IGFBP1","IGFBP2","IGFBP3","IGFBP4","IGFBP5","IGFBP6","IGFBP7","IL10","IL13","IL15","IL18","IL1A","IL1B","IL2","IL32","IL6","IL6ST","IL7","INHA","IQGAP2","ITGA2","ITPKA","JUN","KITLG","LCP1","MIF","MMP1","MMP10","MMP12","MMP13","MMP14","MMP2","MMP3","MMP9","NAP1L4","NRG1","PAPPA","PECAM1","PGF","PIGF","PLAT","PLAU","PLAUR","PTBP1","PTGER2","PTGES","RPS6KA5","SCAMP4","SELPLG","SEMA3F","SERPINB4","SERPINE1","SERPINE2","SPP1","SPX","TIMP2","TNF","TNFRSF10C","TNFRSF11B","TNFRSF1A","TNFRSF1B","TUBGCP2","VEGFA","VEGFC","VGF","WNT16","WNT2"))
names(fgsea_SENMAYO) <- "SenMayo"
pathways.hallmark <- fgsea_SENMAYO

DEG_comp <- res_Short_scr
DEG_comp <- res_ShortvsLong_2
DEG_comp <- res_ShortvsLong_4

DEG_comp <- res_Short_scr_2_paired
DEG_comp <- res_Long_scr_2_paired
DEG_comp <- res_Short_scr_4_paired
DEG_comp <- res_Long_scr_4_paired

#iterate through each comparison for each pathways and save data below.

df <- data.frame(DEG_comp)
df$row <- rownames(DEG_comp)
df$Gene.name<- apply(df,1,function(x) {unlist(strsplit(x[7],"\\|"))[2]})
df2 <- df %>% 
  dplyr::select(Gene.name, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(Gene.name) %>% 
  dplyr::summarize(stat=mean(stat))

ranks <- tibble::deframe(df2)

fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks)

# Tidy the results:
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)

fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()



plotEnrichment(pathways.hallmark[["SenMayo"]],stats=ranks)

#save results after running code above
res_Short_scr_HALLMARK <-data.frame(dataplot)
res_Short_scr_GO <-data.frame(dataplot)
res_Short_scr_SENMAYO <- fgseaResTidy

res_ShortvsLong_2_HALLMARK <-data.frame(dataplot )
res_ShortvsLong_2_GO <-data.frame(dataplot )
res_ShortvsLong_2_SENMAYO <- fgseaResTidy

res_ShortvsLong_4_HALLMARK <-data.frame(dataplot )
res_ShortvsLong_4_GO <-data.frame(dataplot )

res_Short_scr_2_paired_HALLMARK <-data.frame(dataplot)
res_Short_scr_2_paired_GO <-data.frame(dataplot )

res_Long_scr_2_paired_HALLMARK <-data.frame(dataplot )
res_Long_scr_2_paired_GO <-data.frame(dataplot )


res_Short_scr_4_paired_HALLMARK <-data.frame(dataplot)
res_Short_scr_4_paired_GO <-data.frame(dataplot )


res_Long_scr_4_paired_HALLMARK <-data.frame(dataplot )
res_Long_scr_4_paired_GO <-data.frame(dataplot )





##Combine GSEA HALLMARK, Take Top 10 Hallmark NES
df.R<-data.frame(unique(c(res_Short_scr_HALLMARK[order(res_Short_scr_HALLMARK $NES),]$pathway,
                          res_ShortvsLong_2_HALLMARK[order(res_ShortvsLong_2_HALLMARK $NES),]$pathway,
                          res_ShortvsLong_4_HALLMARK[order(res_ShortvsLong_4_HALLMARK $NES),]$pathway,
                          res_Short_scr_2_paired_HALLMARK[order(res_Short_scr_2_paired_HALLMARK $NES),]$pathway,
                          res_Long_scr_2_paired_HALLMARK[order(res_Short_scr_2_paired_HALLMARK $NES),]$pathway,
                          res_Short_scr_4_paired_HALLMARK[order(res_Short_scr_4_paired_HALLMARK $NES),]$pathway,
                          res_Long_scr_4_paired_HALLMARK[order(res_Short_scr_4_paired_HALLMARK $NES),]$pathway)))
colnames(df.R) <-"pathway"
df.R<- merge(df.R,res_Short_scr_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_ShortvsLong_2_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_ShortvsLong_4_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_Short_scr_2_paired_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_Short_scr_4_paired_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_Long_scr_2_paired_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_Long_scr_4_paired_HALLMARK[,c("pathway","NES")],by="pathway", all.x = T)
colnames(df.R)<-c("pathway","All_ShortSc_LongSc","All_Short2_Long2","All_Short4_Long4","All_ShortSc_Short2","All_ShortSc_Short4","All_LongSc_Long2","All_LongSc_Long4")
df.R[is.na(df.R)]<-0
df.R$pathway_simple <- apply(df.R,1,function(x){unlist(strsplit(x[1],"MARK_"))[2]})


df.R2<-data.frame(unique(c(res_Short_scr_HALLMARK[order(res_Short_scr_HALLMARK $padj),]$pathway,
                           res_ShortvsLong_2_HALLMARK[order(res_ShortvsLong_2_HALLMARK $padj),]$pathway,
                           res_ShortvsLong_4_HALLMARK[order(res_ShortvsLong_4_HALLMARK $padj),]$pathway,
                           res_Short_scr_2_paired_HALLMARK[order(res_Short_scr_2_paired_HALLMARK $padj),]$pathway,
                           res_Long_scr_2_paired_HALLMARK[order(res_Short_scr_2_paired_HALLMARK $padj),]$pathway,
                           res_Short_scr_4_paired_HALLMARK[order(res_Short_scr_4_paired_HALLMARK $padj),]$pathway,
                           res_Long_scr_4_paired_HALLMARK[order(res_Short_scr_4_paired_HALLMARK $padj),]$pathway)))
colnames(df.R2) <-"pathway"
df.R2<- merge(df.R2,res_Short_scr_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_ShortvsLong_2_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_ShortvsLong_4_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_Short_scr_2_paired_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_Short_scr_4_paired_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_Long_scr_2_paired_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_Long_scr_4_paired_HALLMARK[,c("pathway","padj")],by="pathway", all.x = T)
colnames(df.R2)<-c("pathway","All_ShortSc_LongSc","All_Short2_Long2","All_Short4_Long4","All_ShortSc_Short2","All_ShortSc_Short4","All_LongSc_Long2","All_LongSc_Long4")
df.R2[is.na(df.R2)]<-1
df.R2$pathway_simple <- apply(df.R2,1,function(x){unlist(strsplit(x[1],"MARK_"))[2]})
df.R2[,c(2:8)]<- -log(df.R2[,c(2:8)])

df.R_M <- reshape2::melt(df.R[,c(9,2:8)])
df.R2_M <- reshape2::melt(df.R2[,c(9,2:8)])
df.R3_M <- merge(df.R_M,df.R2_M, by= c("pathway_simple","variable"))

df.R3_M$variable <- factor(df.R3_M$variable, levels = c("All_ShortSc_LongSc","All_Short2_Long2","All_Short4_Long4","All_ShortSc_Short2","All_ShortSc_Short4","All_LongSc_Long2","All_LongSc_Long4"))  
ggplot(df.R3_M, aes(x=variable, y = pathway_simple, color = value.x, size = value.y)) + 
  geom_point() + 
  scale_color_gradientn(colours =rev(c(brewer.pal(n = 5, name = "PuOr")))) + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  ylab('') +
  theme(axis.ticks = element_blank()) 

ggplot(df.R3_M, aes(x=variable, y = pathway_simple, color = value.x, size = value.y)) + 
  geom_point() + 
  scale_color_gradient2(low = "#E2B8D6", high = "#5f5e60") + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  ylab('') +
  theme(axis.ticks = element_blank()) 



##Combine GSEA GO NES
df.R<-data.frame(unique(c(res_Short_scr_GO[order(res_Short_scr_GO $NES),]$pathway,
                          res_ShortvsLong_2_GO[order(res_ShortvsLong_2_GO $NES),]$pathway,
                          res_ShortvsLong_4_GO[order(res_ShortvsLong_4_GO $NES),]$pathway,
                          res_Short_scr_2_paired_GO[order(res_Short_scr_2_paired_GO $NES),]$pathway,
                          res_Long_scr_2_paired_GO[order(res_Short_scr_2_paired_GO $NES),]$pathway,
                          res_Short_scr_4_paired_GO[order(res_Short_scr_4_paired_GO $NES),]$pathway,
                          res_Long_scr_4_paired_GO[order(res_Short_scr_4_paired_GO $NES),]$pathway)))
colnames(df.R) <-"pathway"
df.R<- merge(df.R,res_Short_scr_GO[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_ShortvsLong_2_GO[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_ShortvsLong_4_GO[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_Short_scr_2_paired_GO[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_Short_scr_4_paired_GO[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_Long_scr_2_paired_GO[,c("pathway","NES")],by="pathway", all.x = T)
df.R<- merge(df.R,res_Long_scr_4_paired_GO[,c("pathway","NES")],by="pathway", all.x = T)
colnames(df.R)<-c("pathway","All_ShortSc_LongSc","All_Short2_Long2","All_Short4_Long4","All_ShortSc_Short2","All_ShortSc_Short4","All_LongSc_Long2","All_LongSc_Long4")
df.R[is.na(df.R)]<-0
df.R$pathway_simple <- apply(df.R,1,function(x){unlist(strsplit(x[1],"GO_"))[2]})


df.R2<-data.frame(unique(c(res_Short_scr_GO[order(res_Short_scr_GO $padj),]$pathway,
                           res_ShortvsLong_2_GO[order(res_ShortvsLong_2_GO $padj),]$pathway,
                           res_ShortvsLong_4_GO[order(res_ShortvsLong_4_GO $padj),]$pathway,
                           res_Short_scr_2_paired_GO[order(res_Short_scr_2_paired_GO $padj),]$pathway,
                           res_Long_scr_2_paired_GO[order(res_Short_scr_2_paired_GO $padj),]$pathway,
                           res_Short_scr_4_paired_GO[order(res_Short_scr_4_paired_GO $padj),]$pathway,
                           res_Long_scr_4_paired_GO[order(res_Short_scr_4_paired_GO $padj),]$pathway)))
colnames(df.R2) <-"pathway"
df.R2<- merge(df.R2,res_Short_scr_GO[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_ShortvsLong_2_GO[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_ShortvsLong_4_GO[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_Short_scr_2_paired_GO[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_Short_scr_4_paired_GO[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_Long_scr_2_paired_GO[,c("pathway","padj")],by="pathway", all.x = T)
df.R2<- merge(df.R2,res_Long_scr_4_paired_GO[,c("pathway","padj")],by="pathway", all.x = T)
colnames(df.R2)<-c("pathway","All_ShortSc_LongSc","All_Short2_Long2","All_Short4_Long4","All_ShortSc_Short2","All_ShortSc_Short4","All_LongSc_Long2","All_LongSc_Long4")
df.R2[is.na(df.R2)]<-1
df.R2$pathway_simple <- apply(df.R2,1,function(x){unlist(strsplit(x[1],"GO_"))[2]})
df.R2[,c(2:8)]<- -log(df.R2[,c(2:8)])

df.R_M <- melt(df.R[,c(9,2:8)])
df.R2_M <- melt(df.R2[,c(9,2:8)])
df.R3_M <- merge(df.R_M,df.R2_M, by= c("pathway_simple","variable"))

df.R3_M$variable <- factor(df.R3_M$variable, levels = c("All_ShortSc_LongSc","All_Short2_Long2","All_Short4_Long4","All_ShortSc_Short2","All_ShortSc_Short4","All_LongSc_Long2","All_LongSc_Long4"))  
ggplot(df.R3_M[grep("EXTRACELLULAR_MATRIX",df.R3_M$pathway_simple),], aes(x=variable, y = pathway_simple, color = value.x, size = value.y)) + 
  geom_point() + 
  scale_color_gradientn(colours =rev(c(brewer.pal(n = 5, name = "PuOr"))),limits=c(-3,3),) + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  ylab('') +
  theme(axis.ticks = element_blank()) #ECM related GO only



ggplot(df.R3_M[grep("EXTRACELLULAR_MATRIX",df.R3_M$pathway_simple),], aes(x=variable, y = pathway_simple, color = value.x, size = value.y)) + 
  geom_point() + 
  scale_color_gradient2(low = "#E2B8D6", high = "#5f5e60", limits=c(-3,3)) + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  ylab('') +
  theme(axis.ticks = element_blank()) #ECM related GO only



save.image("DESeq2_MDS_CD34_bulkRNAseq_GSEA_forPub.Rdata")






#################### Epithelial to Mesenchymal Transistion analysis ##########################
EMT_COMPLEX_GO <- unique(c(
  unlist(res_ShortvsLong_2_HALLMARK[res_ShortvsLong_2_HALLMARK$pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",]$leadingEdge)))

MDS_EMT <-cdds[which(cdds$Gene.name %in% EMT_COMPLEX_GO ==T),]
MDS_EMT <- MDS_EMT[order(MDS_EMT$Gene.name),]
MDS_EMT <-MDS_EMT[,c("GeneID","Gene.name",as.character(meta$Sample))]
MDS_EMT_p <-MDS_EMT[,c("GeneID","Gene.name",as.character(meta_paired$Sample))]
MDS_EMT_scr <-MDS_EMT[,c("GeneID","Gene.name",as.character(meta[meta$Time2 == "Scr",]$Sample))]

library("circlize") ## For color options
name <- MDS_EMT[,c("Gene.name")]
df.OG2 <- data.matrix(log2(MDS_EMT[,c(3:52)]+0.01))
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(MDS_EMT[,c(3:52)])
pheatmap(df.OG2,color = colorRampPalette(c("#f0f0f0",rev(brewer.pal(n = 11, name = "Spectral"))))(100),breaks = c(seq(0,15,by=0.15)),show_rownames = T,cluster_rows =F, cluster_cols=F)

## Generate Z-score for TPM from all samples (USE THIS)###
library(matrixStats)
MDS_EMT1<-log2(MDS_EMT[,c(3:52)]+0.01)
MDS_EMT1_Zscore<- (MDS_EMT1-rowMeans(MDS_EMT1))/(rowSds(as.matrix(MDS_EMT1)))[row(MDS_EMT1)]
name <- MDS_EMT[,c("Gene.name")]
df.OG2 <- data.matrix(MDS_EMT1_Zscore)
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(MDS_EMT[,c(3:52)])
#df.OG2<-df.OG2[,c(as.character(meta[meta$Time != "Screening" & meta$Time2 == "PC2",]$Sample))] #for meta order, ONLY PC2
pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = T,cluster_rows =F, cluster_cols=F)
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),clustering_method="ward.D2",show_rownames = T,cluster_rows =T, cluster_cols=F)


name <- MDS_EMT[,c("Gene.name")]
df.OG2 <- data.matrix(log2(MDS_EMT_scr[,c(3:31)]+0.01))
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(MDS_EMT_scr[,c(3:31)])
pheatmap(df.OG2,color = colorRampPalette(c("#f0f0f0",rev(brewer.pal(n = 11, name = "Spectral"))))(100),breaks = c(seq(0,15,by=0.15)),show_rownames = T,cluster_rows =F, cluster_cols=F)

MDS_EMT1<-log2(MDS_EMT_scr[,c(3:31)]+0.01)
MDS_EMT1_Zscore<- (MDS_EMT1-rowMeans(MDS_EMT1))/(rowSds(as.matrix(MDS_EMT1)))[row(MDS_EMT1)]
name <- MDS_EMT[,c("Gene.name")]
df.OG2 <- data.matrix(MDS_EMT1_Zscore)
row.names(df.OG2) <- name
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),clustering_method="ward.D2",show_rownames = T,cluster_rows =T, cluster_cols=F)


name <- MDS_EMT[,c("Gene.name")]
df.OG2 <- data.matrix(log2(MDS_EMT_scr[,c(3:31)]+0.01))
row.names(df.OG2) <- name
df.OG2a <- data.matrix(log2(MDS_EMT[,c(as.character(meta[meta$Time2 == "PC2",]$Sample))]+0.01))
row.names(df.OG2a) <- name
df.OG2b <- data.matrix(log2(MDS_EMT[,c(as.character(meta[meta$Time2 == "PC4",]$Sample))]+0.01))
row.names(df.OG2b) <- name
df.OG3<-cbind(df.OG2,df.OG2a,df.OG2b) #combine screening and PC2 together into one zscore
pheatmap(df.OG3,color = colorRampPalette(c("#f0f0f0",rev(brewer.pal(n = 11, name = "Spectral"))))(100),breaks = c(seq(0,15,by=0.15)),show_rownames = T,cluster_rows =F, cluster_cols=F)
df.OG3_Zscore<- (df.OG3-rowMeans(df.OG3))/(rowSds(as.matrix(df.OG3)))[row(df.OG3)]
d<-pheatmap(df.OG3_Zscore,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),clustering_method="ward.D2",show_rownames = T,cluster_rows =F, cluster_cols=F)




####EMT Fold Change analysis###
write.table(MDS_EMT_p[,c(3:37)]+0.01,"MDS_normCount_EMT_for_FC_analysis.txt", sep = "\t",col.names =F, row.names=F, quote=F)

#make fold change: awk -v OFS="\t" '{print $1/$1,$2/$1,$3/$3,$4/$3,$5/$5,$6/$5,$7/$5,$8/$8,$9/$8,$10/$8,$11/$11,$12/$11, $13/$13,$14/$13,$15/$13, $16/$16, $17/$16,$18/$16, $19/$19,$20/$19,$21/$21,$22/$21,$23/$21, $24/$24, $25/$24,$26/$24, $27/$27,$28/$27, $29/$27, $30/$30, $31/$30, $32/$32, $33/$32,$34/$34,$35/$34}' MDS_normCount_EMT_for_FC_analysis.txt > MDS_FC_EMT_for_FC_analysis.txt

EMT_FC <- read.table("MDS_FC_EMT_for_FC_analysis.txt", header = F, sep ="\t", stringsAsFactors = F)
colnames(EMT_FC) <- colnames(MDS_EMT_p[,3:37])
EMT_FC <- cbind(MDS_EMT_p[,1:2],EMT_FC)
EMT_FC1 <- EMT_FC[,c(colnames(EMT_FC)[1:2],as.character(meta_paired$Sample))]

name <- MDS_EMT[,c("Gene.name")]
df.OG2 <- data.matrix(log2(EMT_FC[,3:37]))
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(EMT_FC[,3:37])
pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PiYG")))(100),breaks = c(seq(-2,2,by=0.04)),clustering_method="ward.D2",show_rownames = T,cluster_rows =F, cluster_cols=F)






#################### Viral mimicry analysis ##########################
VMG <- read.table("Viral_mimicry_genelist.txt", sep = "\t", header = T, stringsAsFactors = F)

colnames(VMG) <- c("Gene.name","NMgene","category")
Ensb_genelist<-read.table("ENSEMBL_gene_transcript_ID_Genename.txt", sep = "\t", header =T, stringsAsFactors = F)
Ensb_genelist_Ref<-Ensb_genelist[Ensb_genelist$RefSeq.match.transcript != "",]
VMG_Ens_Ref <- merge(VMG, Ensb_genelist_Ref[,c(4,5)], by = "Gene.name")## ONLY REF seq Transcripts
colnames(VMG_Ens_Ref)[4] <- "GeneID"
VMG_Ens <- merge(VMG, Ensb_genelist[,c(4,5)], by = "Gene.name")## ONLY REF seq Transcripts
colnames(VMG_Ens)[4] <- "GeneID"

MDS_VM <-cdds[which(cdds$Gene.name %in% VMG_Ens$Gene.name==T),]
MDS_VM <- MDS_VM[order(MDS_VM$Gene.name),]
MDS_VM <-MDS_VM[,c("GeneID","Gene.name",as.character(meta$Sample))]
MDS_VM_p <-MDS_VM[,c("GeneID","Gene.name",as.character(meta_paired$Sample))]
MDS_VM_scr <-MDS_VM[,c("GeneID","Gene.name",as.character(meta[meta$Time2 == "Scr",]$Sample))]

library("circlize") ## For color options
name <- MDS_VM[,c("Gene.name")]
df.OG2 <- data.matrix(log2(MDS_VM[,c(3:52)]+0.01))
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(MDS_VM[,c(3:52)])
pheatmap(df.OG2,color = colorRampPalette(c("#f0f0f0",rev(brewer.pal(n = 11, name = "Spectral"))))(100),breaks = c(seq(0,15,by=0.15)),show_rownames = T,cluster_rows =F, cluster_cols=F)

## Generate Z-score for TPM from all samples (USE THIS)###
library(matrixStats)
MDS_VM1<-log2(MDS_VM[,c(3:52)]+0.01)
MDS_VM1_Zscore<- (MDS_VM1-rowMeans(MDS_VM1))/(rowSds(as.matrix(MDS_VM1)))[row(MDS_VM1)]
name <- MDS_VM[,c("Gene.name")]
df.OG2 <- data.matrix(MDS_VM1_Zscore)
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(MDS_VM[,c(3:52)])
#df.OG2<-df.OG2[,c(as.character(meta[meta$Time != "Screening" & meta$Time2 == "PC2"& meta$Time2 == "PC4",]$Sample))] #for meta order, ONLY PC2
pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = T,cluster_rows =F, cluster_cols=F)
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),clustering_method="ward.D2",show_rownames = T,cluster_rows =T, cluster_cols=F)


name <- MDS_VM[,c("Gene.name")]
df.OG2 <- data.matrix(log2(MDS_VM_scr[,c(3:31)]+0.01))
row.names(df.OG2) <- name
df.OG2a <- data.matrix(log2(MDS_VM[,c(as.character(meta[meta$Time2 == "PC2",]$Sample))]+0.01))
row.names(df.OG2a) <- name
df.OG2b <- data.matrix(log2(MDS_VM[,c(as.character(meta[meta$Time2 == "PC4",]$Sample))]+0.01))
row.names(df.OG2b) <- name
df.OG3<-cbind(df.OG2,df.OG2a,df.OG2b) #combine screening and PC2 together into one zscore
pheatmap(df.OG3,color = colorRampPalette(c("#f0f0f0",rev(brewer.pal(n = 11, name = "Spectral"))))(100),breaks = c(seq(0,15,by=0.15)),show_rownames = T,cluster_rows =F, cluster_cols=F)
df.OG3_Zscore<- (df.OG3-rowMeans(df.OG3))/(rowSds(as.matrix(df.OG3)))[row(df.OG3)]
d<-pheatmap(df.OG3_Zscore,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),clustering_method="ward.D2",show_rownames = T,cluster_rows =F, cluster_cols=F)




####VM Fold Change analysis###
write.table(MDS_VM_p[,c(3:37)]+0.01,"MDS_normCount_VM_for_FC_analysis.txt", sep = "\t",col.names =F, row.names=F, quote=F)

#make fold change: awk -v OFS="\t" '{print $1/$1,$2/$1,$3/$3,$4/$3,$5/$5,$6/$5,$7/$5,$8/$8,$9/$8,$10/$8,$11/$11,$12/$11, $13/$13,$14/$13,$15/$13, $16/$16, $17/$16,$18/$16, $19/$19,$20/$19,$21/$21,$22/$21,$23/$21, $24/$24, $25/$24,$26/$24, $27/$27,$28/$27, $29/$27, $30/$30, $31/$30, $32/$32, $33/$32,$34/$34,$35/$34}' MDS_normCount_VM_for_FC_analysis.txt > MDS_FC_VM_for_FC_analysis.txt

VM_FC <- read.table("MDS_FC_VM_for_FC_analysis.txt", header = F, sep ="\t", stringsAsFactors = F)
colnames(VM_FC) <- colnames(MDS_VM_p[,3:37])
VM_FC <- cbind(MDS_VM_p[,1:2],VM_FC)
VM_FC1 <- VM_FC[,c(colnames(VM_FC)[1:2],as.character(meta_paired$Sample))]

name <- MDS_VM[,c("Gene.name")]
df.OG2 <- data.matrix(log2(VM_FC[,3:37]))
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(VM_FC[,3:37])
pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PiYG")))(100),breaks = c(seq(-2,2,by=0.04)),clustering_method="ward.D2",show_rownames = T,cluster_rows =F, cluster_cols=F)






#################### IFNY analysis ##########################
IFNy <- unique(c(
  unlist(res_Long_scr_2_paired_HALLMARK[res_Long_scr_2_paired_HALLMARK$pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE",]$leadingEdge),
  unlist(res_Long_scr_4_paired_HALLMARK[res_Long_scr_4_paired_HALLMARK$pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE",]$leadingEdge)))


MDS_IFNy <-cdds[which(cdds$Gene.name %in% IFNy==T),]
MDS_IFNy <- MDS_IFNy[order(MDS_IFNy$Gene.name),]
MDS_IFNy <-MDS_IFNy[,c("GeneID","Gene.name",as.character(meta$Sample))]
MDS_IFNy_p <-MDS_IFNy[,c("GeneID","Gene.name",as.character(meta_paired$Sample))]
MDS_IFNy_scr <-MDS_IFNy[,c("GeneID","Gene.name",as.character(meta[meta$Time2 == "Scr",]$Sample))]

library("circlize") ## For color options
name <- MDS_IFNy[,c("Gene.name")]
df.OG2 <- data.matrix(log2(MDS_IFNy[,c(3:52)]+0.01))
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(MDS_IFNy[,c(3:52)])
pheatmap(df.OG2,color = colorRampPalette(c("#f0f0f0",rev(brewer.pal(n = 11, name = "Spectral"))))(100),breaks = c(seq(0,15,by=0.15)),show_rownames = T,cluster_rows =F, cluster_cols=F)

## Generate Z-score for TPM from all samples (USE THIS)###
library(matrixStats)
MDS_IFNy1<-log2(MDS_IFNy[,c(3:52)]+0.01)
MDS_IFNy1_Zscore<- (MDS_IFNy1-rowMeans(MDS_IFNy1))/(rowSds(as.matrix(MDS_IFNy1)))[row(MDS_IFNy1)]
name <- MDS_IFNy[,c("Gene.name")]
df.OG2 <- data.matrix(MDS_IFNy1_Zscore)
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(MDS_IFNy[,c(3:52)])
#df.OG2<-df.OG2[,c(as.character(meta[meta$Time != "Screening" & meta$Time2 == "PC2",]$Sample))] #for meta order, ONLY PC2
pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = T,cluster_rows =F, cluster_cols=F)
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),clustering_method="ward.D",show_rownames = T,cluster_rows =T, cluster_cols=F)

name <- MDS_IFNy[,c("Gene.name")]
df.OG2 <- data.matrix(log2(MDS_IFNy_scr[,c(3:31)]+0.01))
row.names(df.OG2) <- name
df.OG2a <- data.matrix(log2(MDS_IFNy[,c(as.character(meta[meta$Time2 == "PC2",]$Sample))]+0.01))
row.names(df.OG2a) <- name
df.OG2b <- data.matrix(log2(MDS_IFNy[,c(as.character(meta[meta$Time2 == "PC4",]$Sample))]+0.01))
row.names(df.OG2b) <- name
df.OG3<-cbind(df.OG2,df.OG2a,df.OG2b) #combine screening and PC2 together into one zscore
pheatmap(df.OG3,color = colorRampPalette(c("#f0f0f0",rev(brewer.pal(n = 11, name = "Spectral"))))(100),breaks = c(seq(0,15,by=0.15)),show_rownames = T,cluster_rows =F, cluster_cols=F)
df.OG3_Zscore<- (df.OG3-rowMeans(df.OG3))/(rowSds(as.matrix(df.OG3)))[row(df.OG3)]
d<-pheatmap(df.OG3_Zscore,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),clustering_method="ward.D2",show_rownames = T,cluster_rows =T, cluster_cols=F)




####IFNy Fold Change analysis###
write.table(MDS_IFNy_p[,c(3:37)]+0.01,"MDS_normCount_IFNy_for_FC_analysis.txt", sep = "\t",col.names =F, row.names=F, quote=F)

#make fold change: awk -v OFS="\t" '{print $1/$1,$2/$1,$3/$3,$4/$3,$5/$5,$6/$5,$7/$5,$8/$8,$9/$8,$10/$8,$11/$11,$12/$11, $13/$13,$14/$13,$15/$13, $16/$16, $17/$16,$18/$16, $19/$19,$20/$19,$21/$21,$22/$21,$23/$21, $24/$24, $25/$24,$26/$24, $27/$27,$28/$27, $29/$27, $30/$30, $31/$30, $32/$32, $33/$32,$34/$34,$35/$34}' MDS_normCount_IFNy_for_FC_analysis.txt > MDS_FC_IFNy_for_FC_analysis.txt

IFNy_FC <- read.table("MDS_FC_IFNy_for_FC_analysis.txt", header = F, sep ="\t", stringsAsFactors = F)
colnames(IFNy_FC) <- colnames(MDS_IFNy_p[,3:37])
IFNy_FC <- cbind(MDS_IFNy_p[,1:2],IFNy_FC)
IFNy_FC1 <- IFNy_FC[,c(colnames(IFNy_FC)[1:2],as.character(meta_paired$Sample))]

name <- MDS_IFNy[,c("Gene.name")]
df.OG2 <- data.matrix(log2(IFNy_FC[,3:37]))
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(IFNy_FC[,3:37])
pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PiYG")))(100),breaks = c(seq(-2,2,by=0.04)),clustering_method="ward.D2",show_rownames = T,cluster_rows =F, cluster_cols=F)



########IMMUNE CHECKPOINT GENE ###############
ICG <- read.table("Immune_checkpoint_Gene_list.txt", sep = "\t", header = F, stringsAsFactors = F)
colnames(ICG) <- c("GeneName","GeneIDRef","GeneType","Gene.name")
ICG<- ICG[ICG$GeneType == "stimulatory" | ICG$GeneType == "inhibitory",] 
ICG <- ICG[order(ICG$GeneType,ICG$GeneName),]
neworder <- unique(ICG$Gene.name)

MDS_ICG <-cdds[which(cdds$Gene.name %in% ICG$Gene.name==T),]
MDS_ICG$Gene.name <-factor(MDS_ICG$Gene.name, levels = neworder)
MDS_ICG <- merge(MDS_ICG,ICG[,c("Gene.name","GeneName")], by = "Gene.name")
MDS_ICG <-MDS_ICG[,c("GeneID","Gene.name",as.character(meta$Sample),"GeneName")]
MDS_ICG_p <-MDS_ICG[,c("GeneID","Gene.name",as.character(meta_paired$Sample),"GeneName")]
MDS_ICG_scr <-MDS_ICG[,c("GeneID","Gene.name",as.character(meta[meta$Time2 == "Scr",]$Sample),"GeneName")] 

MDS_ICG <- MDS_ICG[order(MDS_ICG$Gene.name),]
MDS_ICG_p <- MDS_ICG_p[order(MDS_ICG_p$Gene.name),]
MDS_ICG_scr <- MDS_ICG_scr[order(MDS_ICG_scr$Gene.name),]


library("circlize") ## For color options
name <- MDS_ICG[,c("GeneName")]
df.OG2 <- data.matrix(log2(MDS_ICG[,c(3:52)]+0.01))
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(MDS_ICG[,c(3:52)])
pheatmap(df.OG2,color = colorRampPalette(c("#f0f0f0",rev(brewer.pal(n = 11, name = "Spectral"))))(100),breaks = c(seq(0,15,by=0.15)),show_rownames = T,cluster_rows =F, cluster_cols=F)

## Generate Z-score for TPM from all samples (USE THIS)###
library(matrixStats)
MDS_ICG1<-log2(MDS_ICG[,c(3:52)]+0.01)
MDS_ICG1_Zscore<- (MDS_ICG1-rowMeans(MDS_ICG1))/(rowSds(as.matrix(MDS_ICG1)))[row(MDS_ICG1)]
name <- MDS_ICG[,c("GeneName")]
df.OG2 <- data.matrix(MDS_ICG1_Zscore)
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(MDS_ICG[,c(3:52)])
#df.OG2<-df.OG2[,c(as.character(meta[meta$Time != "Screening" & meta$Time2 == "PC2",]$Sample))] #for meta order, ONLY PC2
pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = T,cluster_rows =F, cluster_cols=F)
d<-pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),clustering_method="ward.D",show_rownames = T,cluster_rows =F, cluster_cols=F)

name <- MDS_ICG[,c("Gene.name")]
df.OG2 <- data.matrix(log2(MDS_ICG_scr[,c(3:31)]+0.01))
row.names(df.OG2) <- name
df.OG2a <- data.matrix(log2(MDS_ICG[,c(as.character(meta[meta$Time2 == "PC2",]$Sample))]+0.01))
row.names(df.OG2a) <- name
df.OG2b <- data.matrix(log2(MDS_ICG[,c(as.character(meta[meta$Time2 == "PC4",]$Sample))]+0.01))
row.names(df.OG2b) <- name
df.OG3<-cbind(df.OG2,df.OG2a,df.OG2b) #combine screening and PC2 together into one zscore
pheatmap(df.OG3,color = colorRampPalette(c("#f0f0f0",rev(brewer.pal(n = 11, name = "Spectral"))))(100),breaks = c(seq(0,15,by=0.15)),show_rownames = T,cluster_rows =F, cluster_cols=F)
df.OG3_Zscore<- (df.OG3-rowMeans(df.OG3))/(rowSds(as.matrix(df.OG3)))[row(df.OG3)]
d<-pheatmap(df.OG3_Zscore,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),breaks = c(seq(-2,2,by=0.04)),clustering_method="ward.D2",show_rownames = T,cluster_rows =F, cluster_cols=F)


### Make boxplot distribution of each gene

MDS_ICG_melt <-  MDS_ICG[c(53,3:52)]

MDS_ICG_melt <-reshape2::melt(MDS_ICG_melt)
MDS_ICG_melt$variable <- as.character(MDS_ICG_melt$variable)
colnames(MDS_ICG_melt) <- c("Gene","Sample","normCounts")
MDS_ICG_melt<-merge(meta,MDS_ICG_melt, by="Sample")

MDS_ICG_melt<-MDS_ICG_melt[order(MDS_ICG_melt$Sample),]
MDS_ICG_melt$SubjectID <- factor(MDS_ICG_melt$SubjectID, levels = c(unique(as.character(meta$SubjectID))))

a<-ggplot(MDS_ICG_melt[MDS_ICG_melt$Time2 == "Scr",], aes(x=Survival, y=normCounts)) +
  geom_boxplot(aes(fill=Survival),position = position_dodge(0.9),)+geom_point(aes(y=normCounts,fill=Survival), size=3,pch=21) +ggtitle("Immune Checkpoint Genes (Screening)")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = c("white","grey","#c5a5f3"))+labs(x = "Sample", y = "normCounts")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
a+ facet_wrap( ~ Gene, scales="free",ncol=10)


#Do Welch T for ICG expression
wlt_total <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(wlt_total) <- c("Gene","pval","Comparison")

for (x in unique(MDS_ICG_melt$Gene)) {
  wlt <- t.test(MDS_ICG_melt[MDS_ICG_melt$Gene == x & MDS_ICG_melt$Survival == "Healthy",]$normCounts, MDS_ICG_melt[MDS_ICG_melt$Gene == x & MDS_ICG_melt$Survival == "Short",]$normCounts, alternative = "two.sided", var.equal = FALSE)
  wltdata <- data.frame(x,wlt$p.value,"Healthy_ShortSurv")
  colnames(wltdata) <- c("Gene","pval","Comparison")
  wlt_total <- rbind(wlt_total,wltdata)}
for (x in unique(MDS_ICG_melt$Gene)) {
  wlt <- t.test(MDS_ICG_melt[MDS_ICG_melt$Gene == x & MDS_ICG_melt$Survival == "Healthy",]$normCounts, MDS_ICG_melt[MDS_ICG_melt$Gene == x & MDS_ICG_melt$Survival == "Long",]$normCounts, alternative = "two.sided", var.equal = FALSE)
  wltdata <- data.frame(x,wlt$p.value,"Healthy_LongSurv")
  colnames(wltdata) <- c("Gene","pval","Comparison")
  wlt_total <- rbind(wlt_total,wltdata)}
for (x in unique(MDS_ICG_melt$Gene)) {
  wlt <- t.test(MDS_ICG_melt[MDS_ICG_melt$Gene == x & MDS_ICG_melt$Survival == "Short",]$normCounts, MDS_ICG_melt[MDS_ICG_melt$Gene == x & MDS_ICG_melt$Survival == "Long",]$normCounts, alternative = "two.sided", var.equal = FALSE)
  wltdata <- data.frame(x,wlt$p.value,"ShortSurv_LongSurv")
  colnames(wltdata) <- c("Gene","pval","Comparison")
  wlt_total <- rbind(wlt_total,wltdata)}

wlt_total_p0.05 <- wlt_total[wlt_total$pval < 0.05,]
wlt_total_p0.05$pval <- round(wlt_total_p0.05$pval,3)


####ICG Fold Change analysis###
write.table(MDS_ICG_p[,c(3:37)]+0.01,"MDS_normCount_ICG_for_FC_analysis.txt", sep = "\t",col.names =F, row.names=F, quote=F)

#make fold change: awk -v OFS="\t" '{print $1/$1,$2/$1,$3/$3,$4/$3,$5/$5,$6/$5,$7/$5,$8/$8,$9/$8,$10/$8,$11/$11,$12/$11, $13/$13,$14/$13,$15/$13, $16/$16, $17/$16,$18/$16, $19/$19,$20/$19,$21/$21,$22/$21,$23/$21, $24/$24, $25/$24,$26/$24, $27/$27,$28/$27, $29/$27, $30/$30, $31/$30, $32/$32, $33/$32,$34/$34,$35/$34}' MDS_normCount_ICG_for_FC_analysis.txt > MDS_FC_ICG_for_FC_analysis.txt

ICG_FC <- read.table("MDS_FC_ICG_for_FC_analysis.txt", header = F, sep ="\t", stringsAsFactors = F)
colnames(ICG_FC) <- colnames(MDS_ICG_p[,3:37])
ICG_FC <- cbind(MDS_ICG_p[,1:2],ICG_FC)
ICG_FC1 <- ICG_FC[,c(colnames(ICG_FC)[1:2],as.character(meta_paired$Sample))]

name <- MDS_ICG[,c("GeneName")]
df.OG2 <- data.matrix(log2(ICG_FC[,3:37]))
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(ICG_FC[,3:37])
pheatmap(df.OG2,color = colorRampPalette(rev(brewer.pal(n = 9, name = "PiYG")))(100),breaks = c(seq(-2,2,by=0.04)),show_rownames = T,cluster_rows =F, cluster_cols=F)


ICG_FC_post <- ICG_FC[,-c(grep("SC", colnames(ICG_FC)))]
melt_ICG_FC_post <- melt(ICG_FC_post[,c(2:23)])
colnames(melt_ICG_FC_post) <- c("Gene","Sample","FC")
melt_ICG_FC_post<-merge(meta,melt_ICG_FC_post, by="Sample")


a<-ggplot(melt_ICG_FC_post , aes(x=SurvivalTime, y=log2(FC))) +
  geom_boxplot(aes(fill=SurvivalTime),position = position_dodge(0.9),)+geom_point(aes(y=log2(FC), fill = SurvivalTime), size=2,pch=21) +ggtitle("Immune Checkpoint Genes - Fold Change")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = c("grey","#65737e","#be75be","#872e87"))+labs(x = "Sample", y = "FoldChange")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
a+ facet_wrap( ~ Gene, scales="free",ncol=10)





######################### TE EXPRESSION ANALYSIS ################################
TE_cpms <- read.table("MDS_TE_expression_all_sample_CPM_wHealthyCD34.bed", sep ="\t", header = T, stringsAsFactors = F)

TE_dic <- read.table("GRCh38_GENCODE_rmsk_TE_dictionary.txt",sep ="\t", header = F, stringsAsFactors = F)
colnames(TE_dic) <- c("TE_Subfamily","Geneid","TE_Family","TE_Class")
TE_dic <- TE_dic[!duplicated(TE_dic),]


TE_cpms <- merge(TE_cpms, TE_dic, by = "Geneid", x.all =T)
TE_cpms <- TE_cpms[TE_cpms$TE_Class == "LTR" | TE_cpms$TE_Class == "LINE" | TE_cpms$TE_Class == "SINE"| TE_cpms$TE_Class == "DNA"| TE_cpms$TE_Class == "Retroposon",] #only LINE, SINE, LTR
TE_cpm <- TE_cpms[,c(2,3,4,5,6,1,57,58,59,7:56)]
TE_cpm <- TE_cpm[,c(colnames(TE_cpm[1:9]),as.character(meta$Sample))]
TE_cpm_scr <- TE_cpm[,c(colnames(TE_cpm[1:9]),as.character(meta[meta$Time2 == "Scr",]$Sample))]


TE_annot <- read.delim("GRCh38_GENCODE_rmsk_TE_GTF_parsed.V37annotate.bed", sep ="\t", header =T, stringsAsFactors = F)
TE_annot <- TE_annot[,c(2,3,4,8)]
TE_annot$Anno2 <- apply(TE_annot, 1, function(x) {unlist(strsplit(x[4], " \\("))[1]})
TE_annot$Start <- TE_annot$Start-1 #make into 0-based
TE_cpm_anno_all<- merge(TE_cpm,TE_annot, by=c("Chr","Start","End")) 
TE_cpm_anno_all <-TE_cpm_anno_all[,c(1:9,61,10:59)]
TE_cpm_anno_scr<- merge(TE_cpm_scr,TE_annot, by=c("Chr","Start","End")) 
TE_cpm_anno_scr <-TE_cpm_anno_scr[,c(1:9,40,10:38)]
TE_cpm_anno_scr <-TE_cpm_anno_scr[,c(colnames(TE_cpm_anno_scr)[1:10],as.character(meta[meta$Time2 == "Scr",]$Sample))]


###########Screening Samples only analysis ########################
filter_1cpm <- which(apply(TE_cpm_anno_scr[,11:39], 1, function (x) (max(x)>=1)))#only MDS samples cutoff

TE_cpm_anno_scr1<- TE_cpm_anno_scr[filter_1cpm,]  


TE_cpm_anno_scr1_intergenic <-TE_cpm_anno_scr1[TE_cpm_anno_scr1$Anno2 == "Intergenic",] 
TE_cpm_anno_scr1_intergenic<-TE_cpm_anno_scr1_intergenic[order(TE_cpm_anno_scr1_intergenic$TE_Class,TE_cpm_anno_scr1_intergenic$TE_Family,TE_cpm_anno_scr1_intergenic$TE_Subfamily),] #sort by TE name

TE_cpm_anno_scr1_intergenic_LTR <- TE_cpm_anno_scr1_intergenic[TE_cpm_anno_scr1_intergenic$TE_Class == "LTR",]
TE_cpm_anno_scr1_intergenic_LINE <- TE_cpm_anno_scr1_intergenic[TE_cpm_anno_scr1_intergenic$TE_Class == "LINE",] 


TE_cpm_anno_scr1_intergenic_LTR_melt <-reshape2::melt(TE_cpm_anno_scr1_intergenic_LTR[,c(11:39)])
TE_cpm_anno_scr1_intergenic_LTR_melt$variable <- as.character(TE_cpm_anno_scr1_intergenic_LTR_melt$variable)
colnames(TE_cpm_anno_scr1_intergenic_LTR_melt) <- c("Sample","TotalCPM")
TE_cpm_anno_scr1_intergenic_LTR_melt<-merge(meta,TE_cpm_anno_scr1_intergenic_LTR_melt, by="Sample")

TE_cpm_anno_scr1_intergenic_LTR_melt<-TE_cpm_anno_scr1_intergenic_LTR_melt[order(TE_cpm_anno_scr1_intergenic_LTR_melt$Sample),]
TE_cpm_anno_scr1_intergenic_LTR_melt$SubjectID <- factor(TE_cpm_anno_scr1_intergenic_LTR_melt$SubjectID, levels = c(unique(as.character(meta$SubjectID))))

means_LTR <- aggregate(TotalCPM ~  Sample, TE_cpm_anno_scr1_intergenic_LTR_melt, mean)

TE_cpm_anno_scr1_intergenic_LINE_melt <-reshape2::melt(TE_cpm_anno_scr1_intergenic_LINE[,c(11:39)])
TE_cpm_anno_scr1_intergenic_LINE_melt$variable <- as.character(TE_cpm_anno_scr1_intergenic_LINE_melt$variable)
colnames(TE_cpm_anno_scr1_intergenic_LINE_melt) <- c("Sample","TotalCPM")
TE_cpm_anno_scr1_intergenic_LINE_melt<-merge(meta,TE_cpm_anno_scr1_intergenic_LINE_melt, by="Sample")

TE_cpm_anno_scr1_intergenic_LINE_melt<-TE_cpm_anno_scr1_intergenic_LINE_melt[order(TE_cpm_anno_scr1_intergenic_LINE_melt$Sample),]
TE_cpm_anno_scr1_intergenic_LINE_melt$SubjectID <- factor(TE_cpm_anno_scr1_intergenic_LINE_melt$SubjectID, levels = c(unique(as.character(meta$SubjectID))))

means_LINE <- aggregate(TotalCPM ~  Sample, TE_cpm_anno_scr1_intergenic_LINE_melt, mean)


means_LTR <-merge(meta,means_LTR, by="Sample")
means_LINE <-merge(meta,means_LINE, by="Sample")


a<-ggplot(means_LTR, aes(x=Survival, y=log2(TotalCPM))) +
  geom_boxplot(aes(fill=Survival),position = position_dodge(0.9),)+geom_point(aes(y=log2(TotalCPM), fill=Survival), size=3,pch=21) +ggtitle("Average Intergenic LTR expression Distribution\nin Screening Samples (log2(CPM))")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = c("white","grey","#c5a5f3","black","grey","#c5a5f3"))+labs(x = "Sample", y = "log2 (Ave CPM)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
b<-ggplot(means_LINE, aes(x=Survival, y=log2(TotalCPM))) +
  geom_boxplot(aes(fill=Survival),position = position_dodge(0.9),)+geom_point(aes(y=log2(TotalCPM),fill=Survival), size=3,pch=21) +ggtitle("Average Intergenic LINE expression Distribution\nin Screening Samples (log2(CPM))")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = c("white","grey","#c5a5f3","black","grey","#c5a5f3"))+labs(x = "Sample", y = "log2 (Ave CPM)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")

multiplot(a,b,cols=1)







#########       TE  PAIRED SAMPLES ANALYSIS#######################
filter_1cpm <- which(apply(TE_cpm_anno_pair[,11:45], 1, function (x) (max(x)>=1)))#only MDS samples cutoff

TE_cpm_anno_pair1<- TE_cpm_anno_pair[filter_1cpm,]  

TE_cpm_anno_pair1_intergenic <-TE_cpm_anno_pair1[TE_cpm_anno_pair1$Anno2 == "Intergenic",]
TE_cpm_anno_pair1_intergenic<-TE_cpm_anno_pair1_intergenic[order(TE_cpm_anno_pair1_intergenic$TE_Class,TE_cpm_anno_pair1_intergenic$TE_Family,TE_cpm_anno_pair1_intergenic$TE_Subfamily),] #sort by TE name

TE_cpm_anno_pair1_intergenic_LTR <- TE_cpm_anno_pair1_intergenic[TE_cpm_anno_pair1_intergenic$TE_Class == "LTR",]
TE_cpm_anno_pair1_intergenic_LINE<- TE_cpm_anno_pair1_intergenic[TE_cpm_anno_pair1_intergenic$TE_Class == "LINE",] 



TE_cpm_anno_pair1_intergenic_LTR_melt <-reshape2::melt(TE_cpm_anno_pair1_intergenic_LTR[,c(11:45)])
TE_cpm_anno_pair1_intergenic_LTR_melt$variable <- as.character(TE_cpm_anno_pair1_intergenic_LTR_melt$variable)
colnames(TE_cpm_anno_pair1_intergenic_LTR_melt) <- c("Sample","TotalCPM")
TE_cpm_anno_pair1_intergenic_LTR_melt<-merge(meta,TE_cpm_anno_pair1_intergenic_LTR_melt, by="Sample")

TE_cpm_anno_pair1_intergenic_LTR_melt<-TE_cpm_anno_pair1_intergenic_LTR_melt[order(TE_cpm_anno_pair1_intergenic_LTR_melt$Sample),]
TE_cpm_anno_pair1_intergenic_LTR_melt$SubjectID <- factor(TE_cpm_anno_pair1_intergenic_LTR_melt$SubjectID, levels = c(unique(as.character(meta$SubjectID))))

means_LTR_paired <- aggregate(TotalCPM ~  Sample, TE_cpm_anno_pair1_intergenic_LTR_melt, mean)

TE_cpm_anno_pair1_intergenic_LINE_melt <-reshape2::melt(TE_cpm_anno_pair1_intergenic_LINE[,c(11:45)])
TE_cpm_anno_pair1_intergenic_LINE_melt$variable <- as.character(TE_cpm_anno_pair1_intergenic_LINE_melt$variable)
colnames(TE_cpm_anno_pair1_intergenic_LINE_melt) <- c("Sample","TotalCPM")
TE_cpm_anno_pair1_intergenic_LINE_melt<-merge(meta,TE_cpm_anno_pair1_intergenic_LINE_melt, by="Sample")

TE_cpm_anno_pair1_intergenic_LINE_melt<-TE_cpm_anno_pair1_intergenic_LINE_melt[order(TE_cpm_anno_pair1_intergenic_LINE_melt$Sample),]
TE_cpm_anno_pair1_intergenic_LINE_melt$SubjectID <- factor(TE_cpm_anno_pair1_intergenic_LINE_melt$SubjectID, levels = c(unique(as.character(meta$SubjectID))))

means_LINE_paired <- aggregate(TotalCPM ~  Sample, TE_cpm_anno_pair1_intergenic_LINE_melt, mean)




#distribution
means_LTR_paired <-merge(meta,means_LTR_paired, by="Sample")
means_LTR_paired <-rbind(means_LTR[1:8,],means_LTR_paired) #add healthy

means_LINE_paired <-merge(meta,means_LINE_paired, by="Sample")
means_LINE_paired <-rbind(means_LINE[1:8,],means_LINE_paired) #add healthy

means_LTR_paired$SubjectID <- factor(means_LTR_paired$SubjectID, levels =c(unique(meta$SubjectID)))

means_LINE_paired$SubjectID <- factor(means_LINE_paired$SubjectID, levels =c(unique(meta$SubjectID)))

a<-ggplot(means_LTR_paired, aes(x=SurvivalTime, y=log2(TotalCPM))) +
  geom_boxplot(aes(fill=Survival),position = position_dodge(0.9),)+geom_point(aes(y=log2(TotalCPM),fill=SubjectID), size=4,pch=21) +ggtitle("Average Intergenic LTR expression Distribution (log2(CPM))")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(values = c("white","grey","#c5a5f3","black","black","black","black","black","black","black","black",mypalette,mypalette2))+scale_fill_manual(values = c("white","grey","#c5a5f3","black","black","black","black","black","black","black","black",mypalette,mypalette2))+labs(x = "Sample", y = "log2 (Ave CPM)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
b<-ggplot(means_LINE_paired, aes(x=SurvivalTime, y=log2(TotalCPM))) +
  geom_boxplot(aes(fill=Survival),position = position_dodge(0.9),)+geom_point(aes(y=log2(TotalCPM),fill=SubjectID), size=4,pch=21) +ggtitle("Average Intergenic LINE expression Distribution (log2(CPM))")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(values = c("white","grey","#c5a5f3","black","black","black","black","black","black","black","black",mypalette,mypalette2))+scale_fill_manual(values = c("white","grey","#c5a5f3","black","black","black","black","black","black","black","black",mypalette,mypalette2))+labs(x = "Sample", y = "log2 (Ave CPM)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))

multiplot(a,b, cols =1)




save.image("MDS_TE_expression_forPub.RData")









