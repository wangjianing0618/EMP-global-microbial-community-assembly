bNTI<-function(comm, dis, samp.group=NA, weighted=c(TRUE,FALSE),grouping=c(FALSE,TRUE),rand=1000,output.bMNTD=c(FALSE,TRUE))
{
# calculate betaNTI based on betaMNTD, need package "picante" #
## Improvement from Stegen's paper: you can calculate within group.
## cite the original reference: Stegen JC, Lin X, Fredrickson JK, Chen X, Kennedy DW, Murray CJ et al. (2013). Quantifying community assembly processes and identifying features that impose them. Isme Journal 7: 2069-2079.##
## cite this version of R script as personal communication from Daliang Ning in University of Oklahoma ##

## grouping: If grouping=TRUE, randomization will perform within group. If group="N", randomization will be across all samples. Default is FALSE.##

library(picante) #load package
samp.name=rownames(comm)
result=data.frame(matrix(NA,length(samp.name),length(samp.name)))
colnames(result)=samp.name
rownames(result)=samp.name
result.mntd=result

if(grouping)
{
# calculate betaMNTD and betaNTI within group #
group=levels(as.factor(samp.group[,2]))
group.n=length(group)
bMNTD=list()
bNTI=list()
pbar <- txtProgressBar(min = 0, max = 20, style = 3)# progress bar for large data
	for(m in 1:group.n)
	{
	samp.namex=samp.group[samp.group[,2]==group[m],1]#choose group
	if(length(samp.namex)<2){bNTI[[m]]=NA;bMNTD[[m]]=NA}else{
	comx=comm[match(samp.namex,rownames(comm)),]#remove others
	comx=comx[,colSums(comx != 0, na.rm = TRUE) > 0]#remove undetected species in this group, this part about remove NA, I added it myself
	spnamex=colnames(comx)
	disx=dis[match(spnamex,rownames(dis)),match(spnamex,colnames(dis))]#remove undetected species
	bMNTD.obs<-comdistnt(comx, disx, abundance.weighted = weighted, exclude.conspecifics = FALSE)# calculate observed betaMNTD.
	bMNTD[[m]]=as.matrix(bMNTD.obs)
	bMNTD.rand=array(dim=c(length(samp.namex),length(samp.namex),rand))
		for(i in 1:rand)
		{
		rand.namex=sample(spnamex)
		#disx.rand=dis[match(rand.namex,spnamex),match(rand.namex,spnamex)]
		disx.rand=disx
		colnames(disx.rand)=rand.namex
		rownames(disx.rand)=rand.namex
		bMNTD.rand[,,i]<-as.matrix(comdistnt(comx, disx.rand, abundance.weighted = weighted, exclude.conspecifics = FALSE))
		setTxtProgressBar(pbar, round(((m-1)*rand+i)*20/(group.n*rand),0))
		}
	bNTI[[m]]=(bMNTD[[m]]-apply(bMNTD.rand,c(1,2),mean))/(apply(bMNTD.rand,c(1,2),sd))
									}
	result[match(rownames(bNTI[[m]]),samp.name),match(colnames(bNTI[[m]]),samp.name)]=bNTI[[m]]
	result.mntd[match(samp.namex,samp.name),match(samp.namex,samp.name)]=bMNTD[[m]]
	}
close(pbar)
}else{
# calculate across all samples #
bMNTD.obs<-as.matrix(comdistnt(comm, dis, abundance.weighted = weighted, exclude.conspecifics = FALSE)) # calculate observed betaMNTD.
spname=colnames(comm)
bMNTD.rand=array(dim=c(length(samp.name),length(samp.name),rand))
pbar <- txtProgressBar(min = 0, max = 20, style = 3)# progress bar for large data
	for(i in 1:rand)
	{
	rand.name=sample(spname)
	dis.rand=dis
	colnames(dis.rand)=rand.name
	rownames(dis.rand)=rand.name
	bMNTD.rand[,,i]<-as.matrix(comdistnt(comm, dis.rand, abundance.weighted = weighted, exclude.conspecifics = FALSE))
	setTxtProgressBar(pbar, round((i*20/rand),0))
	}
close(pbar)
bNTI=(bMNTD.obs-apply(bMNTD.rand,c(1,2),mean))/(apply(bMNTD.rand,c(1,2),sd))
result=bNTI
result.mntd=bMNTD.obs
}
if(output.bMNTD)
{
output=list(betaNTI=result,betaMNTD=result.mntd)
}else{
output=result}
output
}
