#set your working directory manually or by this line
setwd("~/Dropbox/Daves_Cnid_Transcriptomes/98-Dave/NV_GO")
system("ls")
rm(list = ls())

#install these at least once, can comment out after first install
#install.packages("ggplot2")
#install.packages("reshape")

#load packages for this session
library(ggplot2); library(reshape)

#load GO data. These can be gotten from column 14 from the .tsv file output of interproscan. Need to specify the -goterms and TSV option. Have a look at my golist.txt files.
#these data are from this paper https://academic.oup.com/icb/article/54/2/276/2797881/Gene-Co-expression-Modules-Underlying-Polymorphic
#allhydrac are the GO terms from the global transcriptome, all others are the modules of co-expressed genes  
#These lines load the data into a dataframe

allNV<-as.data.frame(table(read.table("nvgo_all")$V1))
OF<-as.data.frame(table(read.table("OF_GO_all")$V1)); temp<-merge(allNV,OF,by=1,all.x=T,all.y=T,sort=T)
AG<-as.data.frame(table(read.table("AG_GO_all")$V1)) ;temp<-merge(temp,AG,by=1,all.x=T,all.y=T,sort=T)
Chang<-as.data.frame(table(read.table("chang_GO_all")$V1)); temp<-merge(temp,Chang,by=1,all.x=T,all.y=T,sort=T)
Zap<-as.data.frame(table(read.table("zap_GO_all")$V1));temp<-merge(temp,Zap,by=1,all.x=T,all.y=T,sort=T)

#which is the merged table of GO terms for each sample
GOtable<-temp[,-1]
rownames(GOtable)<-as.character(temp[,1])
colnames(GOtable)<-c("all_NV","OF","AG","Chang","Zappata")

GOtable[is.na(GOtable)]<-0
#GOtable
numsets <-ncol(GOtable)

#This is a function that will annotate with ReviGO. This function is called later and before it will work the .csv files for each GO category (BP, MF, CC) need to be in th
#working directory

makelookup<-function(revigocsv){	
	majorGO<-c(); minorGO<-c()
		for (row in 1:nrow(revigocsv)){
			stat<-revigocsv $plot_X[row]; 
			minorGO[row]<-rownames(revigocsv)[row]
			if (stat!="null"){
				majorGO[row]<-rownames(revigocsv)[row]
				}else{
					majorGO[row]<-majorGO[end(majorGO)[1]]		
					}
			}
			goterms<-cbind(majorGO, minorGO)
			return(as.data.frame(goterms))
}

consolidateCountsbyTerm<-function(GOtablecounts,lookup){
	 sums<-c();
for (m in 1:length(levels(lookup$majorGO)))
{
   numsets<-ncol(GOtablecounts)
	goterm<-levels(lookup$majorGO)[m]
	sum<-c(goterm,colSums(GOtablecounts[GOtablecounts$majorGO==goterm,2: numsets]))
	sums<-rbind(sums, sum)
}
			df<-as.data.frame(sums,row.names=1, stringsAsFactors=F)
	for (i in 1:length(df)) {attr(df[[i]], "names") <- NULL}
	for (c in 2:ncol(df)){df[,c]<-as.numeric(df[,c])}

	return(df)

}


## runs the fisher.test
runFisherTest<-function(GOtab,numsets)
{df <-as.data.frame(GOtab[,1:numsets+1],row.names= GOtab[,1])
fishermatrix<-matrix(NA,nrow=nrow(df),ncol=numsets-1,byrow=T)
for (m in 1:numsets) {
	for (r in 1:nrow(df))
{
geneCounts<-df[r,1:numsets]
totCounts<-colSums(df[,1:numsets])
nonGeneCounts<-totCounts-geneCounts
HydracGOCount<-df[r,1]
HydracOtherGOCount<-totCounts[1]-HydracGOCount
moduleGOCount<-df[r,m]
moduleOtherGOCount<-totCounts[m]-moduleGOCount

pval<-fisher.test(matrix(c(HydracGOCount, HydracOtherGOCount, moduleGOCount, moduleOtherGOCount),nrow = 2))$p.value
fishermatrix[r,m-1]<-pval
}
}
dimnames(fishermatrix)<-list(rownames(df),colnames(df)[2:numsets])
final<-cbind(as.data.frame(fishermatrix),as.character(GOtab $annotation))
return(final)
}

#un-comment the next line if you want to output to pdf
#pdf("GObarplots.pdf",width=11,height=8)

#molfunctiny.csv needs to be in the working dir
molfun<-read.csv("mf_go.csv",header=T, row.names=1)
lookup <-makelookup(molfun)
temp<-GOtable[which(rownames(GOtable) %in% rownames(molfun)),]
majorGO<-as.character(lookup$majorGO)
temp2<-cbind(majorGO,temp[as.character(lookup$minorGO),])
molfunGO <-consolidateCountsbyTerm(temp2,lookup); colnames(molfunGO)[1]<-"terms"
annotation<-molfun[rownames(molfun) %in% molfunGO $terms,"description"]
molfunGO <-cbind(molfunGO, annotation)
molfunGO <-molfunGO[order(molfunGO$annotation),]
#filter out low-gene-count GO terms
#rowSums(molfunGO[,3: numsets+1])>=4
molfunGO <-molfunGO[which(rowSums(molfunGO[,3: numsets+1])>=4),]
molfunGO<-droplevels(molfunGO)

melted<-melt(molfunGO,id=c("terms","annotation"))
#absolute counts
#ggplot(data=melted[melted$variable!="all",], aes(x=variable,y=value,fill=annotation)) + geom_bar(stat="identity")

#frequency
#ggplot(melted) + geom_bar(aes(factor(variable),value, fill=annotation), position='fill',stat="identity") +
#labs(x="Module",y="Proportion of genes represented by GO category", fill="GO: Molecular Function") + theme(legend.position='none') # no legend
#plot it again, to capture legend:
#ggplot(melted) + geom_bar(aes(factor(variable),value, fill=annotation), position='fill',stat="identity") +
#labs(x="Module",y="Proportion of genes represented by GO category", fill="GO: Molecular Function") 

#fisher test
pv<-runFisherTest (molfunGO,numsets)
p.adjusted<-matrix(data=p.adjust(unlist(c(pv[,1:numsets-1]))), nrow=nrow(pv),ncol=ncol(pv)-1,byrow=F)
dimnames(p.adjusted)<-list(rownames(pv),colnames(pv)[1:numsets-1])
annotation<-as.character(molfunGO $annotation)
ptable<-cbind.data.frame(p.adjusted, annotation)
keep<-c()
for (r in 1:nrow(ptable)){
	if (any(ptable[r,1:4]<0.001))
	{keep<-c(keep,r)
		}
	
}
ptable<-ptable[keep,]
molfunGO_sig<-molfunGO[which(molfunGO $annotation %in% ptable$annotation),]

melted<-melt(molfunGO_sig,id=c("terms","annotation"))

#frequency
#ggplot(melted) + geom_bar(aes(factor(variable),value, fill=annotation),   position='fill',stat="identity") + labs(x="Module",y="Proportion of genes #represented by GO category (signif)", fill="GO: Molecular Function") 

ggplot(melted) + geom_bar(aes(factor(variable),value, fill=annotation), position='fill',stat="identity") +
labs(x="Module",y="Proportion of genes represented by GO category (signif)", fill="GO: Molecular Function") + theme(legend.position='none') # no legend

#make relative abundance
molfunGO_ra <-t(molfunGO[,2:6])/colSums(molfunGO[,2:6])
rels<-cbind(molfunGO[,c(1,7)], t(molfunGO_ra))
# look at signif terms
rels[rels$terms %in% molfunGO_sig$terms,]
molfunGO[molfunGO$terms %in% rownames(ptable),]
ptable



####cellcomptiny.csv needs to be in the working dir
cellcomp<-read.csv("cc_go.csv",header=T,row.names=1)
lookup <-makelookup(cellcomp)
temp<-GOtable[which(rownames(GOtable) %in% rownames(cellcomp)),]
majorGO<-as.character(lookup$majorGO)
temp2<-cbind(majorGO,temp[as.character(lookup$minorGO),])
ccGO <-consolidateCountsbyTerm(temp2,lookup); colnames(ccGO)[1]<-"terms"
annotation<-cellcomp[rownames(cellcomp) %in% ccGO $terms,"description"]
ccGO <-cbind(ccGO, annotation)
ccGO <-ccGO[order(ccGO $annotation),]
#filter out low-gene-count GO terms
#rowSums(ccGO[,3: numsets+1])>=4
ccGO <-ccGO[which(rowSums(ccGO[,3: numsets+1])>=4),]
ccGO <-droplevels(ccGO)

melted<-melt(ccGO,id=c("terms","annotation"))

#absolute counts
#ggplot(data=melted[melted$variable!="all",], aes(x=variable,y=value,fill=annotation)) + geom_bar(stat="identity")
#+ theme(legend.position='none') # no legend
#frequency
#ggplot(melted) + geom_bar(aes(factor(variable),value, fill=annotation),   position='fill',stat="identity") + labs(x="Module",y="Proportion of genes represented by GO category", fill="GO: Cellular Component") 
## fisher.test
pv<-runFisherTest(ccGO,numsets)
p.adjusted<-matrix(data=p.adjust(unlist(c(pv[,1:numsets-1]))), nrow=nrow(pv),ncol=ncol(pv)-1,byrow=F)
dimnames(p.adjusted)<-list(rownames(pv),colnames(pv)[1:numsets-1])
annotation<-as.character(ccGO $annotation)
ptable<-cbind.data.frame(p.adjusted, annotation)
keep<-c()
for (r in 1:nrow(ptable)){
	if (any(ptable[r,1:4]<0.001))
	{keep<-c(keep,r)
		}
	
}
ptable<-ptable[keep,]
ccGO_sig<-ccGO[which(ccGO $annotation %in% ptable$annotation),]

melted<-melt(ccGO_sig,id=c("terms","annotation"))

#frequency
#ggplot(melted) + geom_bar(aes(factor(variable),value, fill=annotation),   position='fill',stat="identity") + labs(x="Module",y="Proportion of genes #represented by GO category (signif)", fill="GO: Cell Component") 

ggplot(melted) + geom_bar(aes(factor(variable),value, fill=annotation), position='fill',stat="identity") +
labs(x="Module",y="Proportion of genes represented by GO category (signif)", fill="GO: Molecular Function") + theme(legend.position='none') # no legend

#make relative abundance
ccGO_ra <-t(ccGO[,2:6])/colSums(ccGO[,2:6])
rels<-cbind(ccGO[,c(1,7)], t(ccGO_ra))
# look at signif terms
rels[rels$terms %in% ccGO_sig$terms,]
ccGO[ccGO$terms %in% rownames(ptable),]
ptable


### bioproctiny.csv needs to be in working dir
bioproc<-read.csv("bp_go.csv",header=T,row.names=1)
#120 rows
lookup<-makelookup(bioproc)

temp<-GOtable[which(rownames(GOtable) %in% rownames(bioproc)),]
majorGO<-as.character(lookup$majorGO)
temp2<-cbind(majorGO,temp[as.character(lookup$minorGO),])
bioprocGO <-consolidateCountsbyTerm(temp2,lookup); colnames(bioprocGO)[1]<-"terms"
annotation<-bioproc[rownames(bioproc) %in% bioprocGO$terms,"description"]
bioprocGO<-cbind(bioprocGO, annotation)
bioprocGO <-bioprocGO[order(bioprocGO $annotation),]
#filter out low-gene-count GO terms
#rowSums(bioprocGO[,3: numsets+1])>=4
bioprocGO<-bioprocGO[which(rowSums(bioprocGO[,3: numsets+1])>=4),]
bioprocGO <-droplevels(bioprocGO)
melted<-melt(bioprocGO,id=c("terms","annotation"))

#absolute counts
#ggplot(data=melted[melted$variable!="all",], aes(x=variable,y=value,fill=annotation)) + geom_bar(stat="identity")

#frequency
#ggplot(melted) + geom_bar(aes(factor(variable),value, fill=annotation),   position='fill',stat="identity") + labs(x="Module",y="Proportion of genes represented by GO category", fill="GO: Biological Process") 




##fisher test
pv<-runFisherTest (bioprocGO,numsets)
p.adjusted<-matrix(data=p.adjust(unlist(c(pv[,1:numsets-1]))), nrow=nrow(pv),ncol=ncol(pv)-1,byrow=F)
dimnames(p.adjusted)<-list(rownames(pv),colnames(pv)[1:numsets-1])
annotation<-as.character(bioprocGO $annotation)
ptable<-cbind.data.frame(p.adjusted, annotation)
keep<-c()
for (r in 1:nrow(ptable)){
	if (any(ptable[r,1:4]<0.001))
	{keep<-c(keep,r)
		}
	
}
ptable<-ptable[keep,]
bioprocGO_sig<-bioprocGO[which(bioprocGO$annotation %in% ptable$annotation),]
nrow(bioprocGO_sig)

melted<-melt(bioprocGO_sig,id=c("terms","annotation"))

ptable<-cbind(rownames(ptable),ptable); all_NV<-rep(1,nrow(ptable));
ptable<- cbind(ptable[,1],all_NV,ptable[,2:6]); colnames(ptable)[1]<-"term"
 pmelt<-melt(ptable, id=c("term","annotation")); 
 signif<-pmelt$value<0.001
 signif[pmelt$value<0.001]<-0
 signif[pmelt$value>0.001]<-1
melt2<-cbind(melted, factor(signif));colnames(melt2)[5]<-"pval"

bioprocGO[bioprocGO$terms %in% rownames(ptable),]

#frequency
#ggplot(melted) + 
#geom_bar(aes(factor(variable),value, fill=annotation),   position='fill',stat="identity") + 
#labs(x="Module",y="Proportion of genes represented by GO category (signif)", fill="GO: Biological Process") 

ggplot(melted) + geom_bar(aes(factor(variable),value, fill=annotation), position='fill',stat="identity") +
labs(x="Module",y="Proportion of genes represented by GO category (signif)", fill="GO: Molecular Function") + theme(legend.position='none') # no legend

#with signif noted 

ggplot(melt2) + 
geom_bar(aes(factor(variable),value, fill=signif),   position='fill',stat="identity") + labs(x="Module",y="Proportion of genes represented by GO #category (signif)", fill="GO: Biological Process") 

#make relative abundance
bioproc_ra <-t(bioprocGO[,2:6])/colSums(bioprocGO[,2:6])
rels<-cbind(bioprocGO[,c(1,7)], t(bioproc_ra))
# look at signif terms
rels[rels$terms %in% bioprocGO_sig$terms,]
bioprocGO[bioprocGO$terms %in% rownames(ptable),]
ptable


#un-comment the next line if you want to output to pdf
#dev.off()
#library(gplots)
#heatmap.2(t(t(bioprocGO_sig[,2:6]))
#heatmap.2(t(t(bioprocGO_sig[,2:6])),scale="row")