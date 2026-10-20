# install "dplyer" package
#install.packages("dplyer")
# Load the dplyr package

#Setup your working directory
setwd("C:\\Users\\ASUS\\Desktop\\New folder (2)\\GSE50006")

library(dplyr)

# Load the dataset
results=read.csv("CLL2.csv", header=TRUE)
head(results)
dim(results)

##Select Upregulated genes
Upregulated_data<- results %>% filter(  adj.P.Val <0.05 & logFC >0.8)
dim(Upregulated_data)
head(Upregulated_data)
UPDE <- dplyr::filter(Upregulated_data, !is.na(Gene.symbol))
tail(UPDE )
dim(UPDE )
write.csv(UPDE ,file="UPDE_CLL_1.csv")

##Select Downregulated genes
Downregulated_data <- results %>% filter(  adj.P.Val <0.05 & logFC < -0.8)
dim(Downregulated_data)
head(Downregulated_data)
DownDE <- dplyr::filter(Downregulated_data, !is.na(Gene.symbol))
tail(DownDE)
dim(DownDE)
write.csv(DownDE,file="DownDE_CLL_1.csv")

#################################
DownDEN <- read.csv("DownDE_CLL_1.csv",header=T)
length(unique(DownDEN$Gene.symbol))
DE_D <- unique(DownDEN$Gene.symbol)

UPDEN <- read.csv("UPDE_CLL_1.csv",header=T)
length(unique(UPDEN$Gene.symbol))
DE_U <- unique(UPDEN$Gene.symbol)

DEG <- c(DE_D,DE_U)
length(DEG)
write.csv(DEG,file="DEG_CLL_2.csv")

