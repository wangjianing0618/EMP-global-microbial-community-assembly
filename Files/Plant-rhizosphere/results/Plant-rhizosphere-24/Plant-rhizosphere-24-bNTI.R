library(picante)
source("Files/raw/bNTI.R")
comm1<-read.table("Files/Plant-rhizosphere/results/Plant-rhizosphere-24/Plant-rhizosphere-24-t.txt", head=T)
colnames(comm1)
colnames(comm1)=sub("\\_","-", colnames(comm1)) 
comm1[is.na(comm1)] <- 0
comm1<-t(comm1)
tree<-read.tree(file="Files/Plant-rhizosphere/results/Plant-rhizosphere-24/Plant-rhizosphere-24.nwk")
samp.group<-read.table("Files/Plant-rhizosphere/results/Plant-rhizosphere-24/Plant-rhizosphere-24.group.txt",header=TRUE) 
dis<-cophenetic(tree)
spname<-rownames(dis)
colnames(dis)=rownames(dis)
comm<-comm1[match(spname,rownames(comm1)),1:ncol(comm1)]
rownames(comm)<-spname
comm<-comm[,match(samp.group[,1],colnames(comm))]
comm<-t(comm)
bNTI.group=bNTI(comm, dis, samp.group=samp.group, weighted=TRUE,grouping=TRUE,rand=1000,output.bMNTD=TRUE)
write.csv(bNTI.group$betaNTI,file=" Files/Plant-rhizosphere/results/Plant-rhizosphere-24/bNTI.Plant-rhizosphere-24.csv")
write.csv(bNTI.group$betaMNTD,file=" Files/Plant-rhizosphere/results/Plant-rhizosphere-24/bMNTD.Plant-rhizosphere-24.csv")
