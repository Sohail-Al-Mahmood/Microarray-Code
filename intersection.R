D1=read.csv("DEG_CLL_1.csv",header = TRUE)
D2=read.csv("DEG_CLL_2.csv", header = TRUE)
#D3=read.csv("CLL & T2D_f  DEG.csv", header = TRUE)

T1=D1$x.1
T2=D2$x.1
#T3=D3$x

Unique=Reduce(intersect,list(T1,T2))
length(Unique)
write.csv(Unique,"Combined CLL.csv")


####################################

rm(list = ls())       # Remove all data

# Set working directory and read data
Data1 <- read.csv(file="GSE29226_series_matrix.csv", row.names = 1)
Data2 <- read.csv(file="GSE29221_series_matrix.csv", row.names = 1)


Data1<-Data1[,order(colnames(Data1),decreasing = T)]
Data2<-Data2[,order(colnames(Data2),decreasing = T)]

log2_transformation = function(dataFrame){
  ex <- dataFrame
  qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) { ex[which(ex <= 0)] <- NaN
  dataFrame = log2(ex)
  }
  return(dataFrame)
}
Data1[Data1 <= 0] <- NA  # Replace non-positive values with NA
Data2[Data2 <= 0] <- NA  # Replace non-positive values with NA

Data1 = log2_transformation(Data1)
Data2 = log2_transformation(Data2)

normalization = function(dataFrame){
  d<-dim(dataFrame)
  new_df<-matrix(0,d[1],d[2])
  med_1<-apply(dataFrame,1,median)
  IQR_1<-apply(dataFrame,1,IQR)
  for(i in 1:d[1]){new_df[i,]<-(as.numeric(dataFrame[i,])-med_1[i])/IQR_1[i]}
  rownames(new_df)<-rownames(dataFrame)
  colnames(new_df)<-colnames(dataFrame)
  return(new_df)
}
Data1 <- na.omit(Data1)      # Remove rows with NA values
Data2 <- na.omit(Data2)      # Remove rows with NA values

Data1 = normalization(Data1)
Data2 = normalization(Data2)



common<-sort(Reduce("intersect",list(rownames(Data1),rownames(Data2))))
Data1_filtered<-Data1[common,]
Data2_filtered<-Data2[common,]


combined<-cbind(Data1_filtered,Data2_filtered)
combined<-combined[,order(colnames(combined), decreasing = T)]
View(combined)
length(combined)
write.csv(combined,"COMBINEDT2D_f.csv")
#colnames(combined)<-c(rep("control",18), rep("case",25))
#View(combined)

library(limma)
library(gclus)
library(Biobase)
library(dplyr)
library(cluster)
Data=read.csv(file="COMBINEDT2D_f.csv", row.names = 1)
ControlCase <- c(rep('control', 24), rep('case', 24))
group <- factor(ControlCase)
design <-model.matrix(~ group)
fit <- lmFit(Data, design)
fit2 <- eBayes(fit)
results <- topTable(fit2, coef=2,  adjust="BH", n = nrow(Data), sort.by = "none")
head(results)
write.csv(results,file="ResultT2D_f.csv")


##Select Upregulated genes

UP_Data <- results %>% filter(adj.P.Val <0.05& logFC >0.9)
dim(UP_Data )
UP_Genes <-row.names(UP_Data);UP_Genes
write.csv(UP_Genes,file="UPT2D_f.csv")


##Select Downregulated genes

Down_Data <- results %>% filter(  adj.P.Val <0.05 & logFC < -0.9)
dim(Down_Data )
Down_Genes <-row.names(Down_Data);Down_Genes
write.csv(Down_Genes,file="DownT2D_f.csv")

DEGs <- c(UP_Genes,Down_Genes)
write.csv(DEGs,file="DEGT2D_f.csv")
dim(DEGs)
####################
Gene <- read.csv("DEGs_GeneSymbol.csv")
DEGs_Final <- unique(Gene$To)
write.csv(DEGs_Final,file="DEGs_Final.csv")