rm(list=ls())
# Set Working Directory
setwd("C:\\Users\\ASUS\\Desktop\\Puspo docking\\Siam_p\\Puspo_docking_result\\NFIC\\NFIC_result")

## Read the all dataset
filelist = list.files(pattern = ".*.log")
n <- length(filelist)
n
#assuming tab separated values with a header    
dock = lapply(filelist, 
              function(x)read.table(x, sep = ",", fill = TRUE ,header = FALSE)
) 
#strsplit(strsplit(dock[[172]][2,1]," ")[[1]][4],"_")[[1]][1]
Drug <- c()
for(i in 1:n){
  Drug[i] <- strsplit(strsplit(dock[[i]][2,1]," ")[[1]][4],"_out")[[1]][1] 
}
Drug

#strsplit(dock[[172]][14,1]," ")[[1]][13]
NFIC<- c()
for(i in 1:n){
  A1 <- strsplit(dock[[i]][14,1]," ")[[1]][13]
  A2 <- strsplit(dock[[i]][14,1]," ")[[1]][12]
  NFIC[i] <- paste(A1, A2)
}
NFIC
result <- cbind(Drug,NFIC)
write.csv(result,file="NFIC_result.csv")
