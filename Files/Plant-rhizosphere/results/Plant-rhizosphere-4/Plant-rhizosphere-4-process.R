library(vegan)
library(parallel)
source("/mnt/sdb/wjn/EMP/0-NULL/Plant-rhizosphere/results/Plant-rhizosphere-4/RC.p.R")
com<-read.table("/mnt/sdb/wjn/EMP/0-NULL/Plant-rhizosphere/results/Plant-rhizosphere-4/Plant-rhizosphere-4-t.txt",header=T,sep="",row.names=1) 
comm<-com[,colSums(com)>0]
BC<-vegdist(comm,method="bray")
write.csv(as.matrix(BC),"/mnt/sdb/wjn/EMP/0-NULL/Plant-rhizosphere/results/Plant-rhizosphere-4/Plant-rhizosphere-4-BC.csv")
rc<-RC.p(comm,method="bray",rand=1000,portion=FALSE,nworker=4,memory.G=50)
write.csv(rc,"/mnt/sdb/wjn/EMP/0-NULL/Plant-rhizosphere/results/Plant-rhizosphere-4/Plant-rhizosphere-4-RC.csv")
bNTI<-read.table("/mnt/sdb/wjn/EMP/0-NULL/Plant-rhizosphere/results/Plant-rhizosphere-4/bNTI.Plant-rhizosphere-4.csv",header=T,sep=",",row.names=1)
colnames(bNTI)=rownames(bNTI)
rcc=rc[match(rownames(bNTI),rownames(rc)),match(colnames(bNTI),colnames(rc))]
bNTI.v<-as.vector(as.dist(bNTI))
rc.v<-as.vector(as.dist(rcc))
id.selectna<-(bNTI.v<=2&bNTI.v>=(-2))
num.pair<-length(bNTI.v)
select.h<-sum(bNTI.v>2)/num.pair
select.l<-sum(bNTI.v<(-2))/num.pair
disper.h<-sum(rc.v[id.selectna]>0.95)/num.pair
disper.l<-sum(rc.v[id.selectna]<(-0.95))/num.pair
drift<-sum(rc.v[id.selectna]<=0.95&rc.v[id.selectna]>=(-0.95))/num.pair
res=data.frame(select.h,select.l,disper.h,disper.l,drift,num.pair)
write.csv(res,"/mnt/sdb/wjn/EMP/0-NULL/Plant-rhizosphere/results/Plant-rhizosphere-4/Plant-rhizosphere-4-Processes.csv")
library(vegan)
library(parallel)
source("Files/raw/RC.p.R")
com<-read.table("Files/Plant-rhizosphere/results/Plant-rhizosphere-4/Plant-rhizosphere-4-t.txt",header=T,sep="",row.names=1) 
comm<-com[,colSums(com)>0]
BC<-vegdist(comm,method="bray")
write.csv(as.matrix(BC)," Files/Plant-rhizosphere/results/Plant-rhizosphere-4/Plant-rhizosphere-4-BC.csv")
rc<-RC.p(comm,method="bray",rand=1000,portion=FALSE,nworker=4,memory.G=50)
write.csv(rc," Files/Plant-rhizosphere/results/Plant-rhizosphere-4/Plant-rhizosphere-4-RC.csv")
bNTI<-read.table(" Files/Plant-rhizosphere/results/Plant-rhizosphere-4/bNTI.Plant-rhizosphere-4.csv",header=T,sep=",",row.names=1)
colnames(bNTI)=rownames(bNTI)
rcc=rc[match(rownames(bNTI),rownames(rc)),match(colnames(bNTI),colnames(rc))]
bNTI.v<-as.vector(as.dist(bNTI))
rc.v<-as.vector(as.dist(rcc))
id.selectna<-(bNTI.v<=2&bNTI.v>=(-2))
num.pair<-length(bNTI.v)
select.h<-sum(bNTI.v>2)/num.pair
select.l<-sum(bNTI.v<(-2))/num.pair
disper.h<-sum(rc.v[id.selectna]>0.95)/num.pair
disper.l<-sum(rc.v[id.selectna]<(-0.95))/num.pair
drift<-sum(rc.v[id.selectna]<=0.95&rc.v[id.selectna]>=(-0.95))/num.pair
res=data.frame(select.h,select.l,disper.h,disper.l,drift,num.pair)
write.csv(res," Files/Plant-rhizosphere/results/Plant-rhizosphere-4/Plant-rhizosphere-4-Processes.csv")
