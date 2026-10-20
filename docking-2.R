rm(list=ls())

#setwd("F:/diseases/Docking")
library(gdata)
library("readxl")
 library(plotrix)
######################################################
y <- as.data.frame(read.csv("Top_Ranked_30_drugs.csv", header=T, row.names=1))
y1 <- as.data.frame(apply(y, 2, as.numeric))  # Convert all variable types to numeric
sapply(y1, class) 
dim(y1)
#row.names(y)
y1=y1[1:30]
png(file = 'F:/diseases/Docking/30newDRUGS12.png', width=3000, height=1500, units='px', res=400, bg='white')
#View(x)

####
par(mar = c(6, 4, 1, 8)) # Set the margin on all sides to 2
cellcol<-matrix(rep("#000",270),nrow=9)
cellcol[y1<(-6.99)]<-color.scale(y1[y1<(-6.99)],c(1,0.8),c(0,0.5),0)
cellcol[y1>-6.99]<-color.scale(y1[y1>-6.99],0,c(0.7,0.8),c(0.5,0))
color2D.matplot(x=round(y1,2),cellcolors=cellcol,xlab="",
                ylab="", axes = FALSE, vcex = 2, vcol = "blue",) #,show.values=TRUE
mtext(text = "RECEPTORS PROTEINS (PROPOSED)",
      side = 2, #side 2 = left
      line = 3, cex=.8)
mtext(text = "TOP RANKED 30 META-DRUG AGENTS OUT OF 207 ",
      side = 1, #side 2 = left
      line = 4.5, cex=.7)

legval<-seq(min(y1),max(y1),length.out=6)
#legval<-c(-13 -11 -9  -10  -9  -7)
legcol<-rep("#000",6)
legcol[legval<(6.99)]<-color.scale(legval[legval<(6.99)],c(1,0.8),c(0,0.8),0)
legcol[legval>-6.99]<-color.scale(legval[legval>-6.99],0,c(0.7,0.8),c(0.5,0))
color.legend(30.9,1,31.1,6,legval, rect.col=legcol,align="rb",gradient="y", cex=0.5)


axis(1, at = seq_len(ncol(y1)) - 0.5,labels = FALSE, tck=-0.02)
text(seq(0.5, ncol(y1), by = 1), par("usr")[3]-0.32, srt = 90, adj = 1, xpd = TRUE, labels = names(y1), cex = .4,col = c("red","red","red","black","black", "black","black","red","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black"))
lablist.y<-as.vector(c("GATA2","FOXC1","TKT","SCARB1","PPA1","GSN","FYN","FKBP1A","ACTN1"))  
axis(2, at=seq(0.5, 9, by=1), labels = FALSE, tck=-0.02)
text(y = seq(0.5, 9, by=1), par("usr")[3]+0.05, labels = lablist.y, srt = 0, pos = 2, xpd = TRUE, col = c("blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue"), cex = .6)

mtext("Binding Affinity",side=4,line=-.3, 
      at=par("usr")[1]+0.17*diff(par("usr")[1:2]),
      cex=0.6)

dev.off()
