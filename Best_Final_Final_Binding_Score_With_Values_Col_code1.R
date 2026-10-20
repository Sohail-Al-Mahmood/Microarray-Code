rm(list=ls())

setwd("C:\\Users\\ASUS\\Desktop\\Puspo docking\\docking_result")
library("gdata")
library("readxl")
library("plotrix")



######################################################
y <- as.data.frame(read.csv("Docking_final_result.csv", header=T, row.names=1))
y1 <- as.data.frame(apply(y, 2, as.numeric))  # Convert all variable types to numeric
sapply(y1, class) 
dim(y1)
#row.names(y)
#View(row.names(y))

png(file ='C:\\Users\\ASUS\\Desktop\\Puspo docking\\docking_result\\Final_Heatmap_35 Drugs.png', width=3000, height=1500, units='px', res=400, bg='white') # pointsize=2 if we don't use vcex in color2D.matplot


#prepare heatmap sides [par(mar = c(bottom, left, top, right))] 
par(mar = c(5, 3.5, 1.5, 3)) 
#  par(mar = c(50, 30, 5, 30)) in case of pointsize=2 
cellcol<-matrix(rep("#000000",315),nrow=9)
cellcol[y1<(-6.99)]<-color.scale(y1[y1<(-6.99)],c(1,0.8),c(0,0.5),0)
cellcol[y1>(-6.99)]<-color.scale(y1[y1>(-6.99)],0,c(0.7,0.8),c(0.5,0))
color2D.matplot(x=round(y1,2),cellcolors=cellcol,xlab="",
                ylab="", axes = FALSE, vcex = 0.3, vcol = "black", show.values=2)
mtext(text = "Top-ranked proposed receptors",
      side = 2.8, #side =2 = left
      line = 2.4, cex=0.8, adj=0.5) #cex=1=default
mtext(text = "Top-ranked 35 meta-drug agents out of 310 drugs",
      side = 1.8, #side =1 = Bottom
      line = 4, cex=0.8, adj=0.5)



legval<-seq(min(y1),max(y1),length.out=6)
#legval<-c(-13 -11 -9  -10  -9  -7)
legcol<-rep("#000",6)
legcol[legval<(6.99)]<-color.scale(legval[legval<(6.99)],c(1,0.8),c(0,0.8),0)
legcol[legval>-6.99]<-color.scale(legval[legval>-6.99],0,c(0.7,0.8),c(0.5,0))
#prepare label in right side
color.legend(36.1,1,36.3,8,legval, rect.col=legcol,align="rb",gradient="y", cex=0.5)



#legval<-seq(min(y1),max(y1),length.out=6)
#legval<-c(-10.2 -9.3 -8.4 -7.5  -6.6  -5.7) ####by=0.9
#legcol<-rep("#000000",6)
#legcol[legval<(-6.99)]<-color.scale(legval[legval<(-6.99)],c(1,0.8),c(0,0.5),0)
#legcol[legval>(-6.99)]<-color.scale(legval[legval>(-6.99)],0,c(0.7,0.8),c(0.5,0))
#color.legend(30.8,2,35.1,8,legval, rect.col=legcol,align="rb",gradient="y", cex=0.5)


#axis(1, at = seq_len(ncol(y1))-0.5,labels = FALSE, tck=-0.02)
#text(seq(0.5, ncol(y1), by = 1), par("usr")[3]-0.32, srt = 45, adj = 1, xpd = TRUE, labels = names(y1),col = c("red","red","red","black","black","black","black","red", "black","black","black","black","black","black","black","black" ,"black" ,"black","black","black","black","black", "black", "black","black","black","black", "black","black","black" ), cex = 0.6) #srt=90 = Drug names straight
#lablist.y<-as.vector(c("GATA2","FOXC1","TKT","SCARB1","PPA1","GSN","FYN","FKBP1A","ACTN1"))    
#axis(2, at=seq(0.5, 9, by=1), labels = FALSE, tck=-0.02)
#text(y=seq(0.5, 9, by=1), par("usr")[3]+0.05, labels = lablist.y, srt = 0, pos = 2, xpd = TRUE, col = c("black","black","black","black","black","black","black","red", "red","red"), cex = 0.5)

axis(1, at = seq_len(ncol(y1)) - 0.5,labels = FALSE, tck=-0.02)
text(seq(0.5, ncol(y1), by = 1), par("usr")[3]-0.32, srt = 90, adj = 1, xpd = TRUE, labels = names(y1), cex = .55,col = c("red","blue","red","blue","red", "black","red","red","red","red","black","blue","black","blue","red","red","red","black","blue","blue","red","black","black","black","black","black","red","black","black","black","red","black","black","red","black"))
lablist.y<-as.vector(c("FOXC1","GATA2","NFIC","BUB1B","KIF11","CCNB2","CCNA2","CCNB1","CDK1"))  
axis(2, at=seq(0.5, 9, by=1), labels = FALSE, tck=-0.02)
text(y = seq(0.5, 9, by=1), par("usr")[3]+0.05, labels = lablist.y, srt = 0, pos = 2, xpd = TRUE, col = c("blue","blue","red","blue","red","red","red","red","red"), cex = .55)




#mtext(text = "Binding Affinity Scores",
      #side = 4, #side 4 = Right
      #line =-0.3, cex=0.6, adj=0.5) #  at=par("usr")[1]+0.18*diff(par("usr")[1:2]),

mtext("Binding Affinity Scores",side=4,line=-.1, 
      at=par("usr")[1]+0.17*diff(par("usr")[1:2]),
      cex=0.6)


dev.off()





















