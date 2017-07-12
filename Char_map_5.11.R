setwd("~/Dropbox/Daves_Cnid_Transcriptomes/98-Dave/character_mapping/2.20_charmap/")
#setwd("~/Desktop/charmap/")
system("ls")

	#install.packages("phytools")
	require(phytools)
	packageVersion("phytools")
	
#library(colorspace); rainbow_hcl(4)  #optional package to pick nice palettes
colshcl<-c("#E495A5" ,"#ABB065" ,"#39BEB1", "#ACA4E2","snow3")
	#par(mfrow=c(1,4))
	
	tr<-ladderize(read.nexus("Cnid_only.nex"))
	tr$tip.label
	tr$tip.label<-unlist(strsplit(tr$tip.label,split="_v1"))[c(T,F)]
	#tree<-root(tr, outgroup="Mnemiopsis_leidyi")
	tree<-tr
	#plot(tree, cex=0.7)
	
	tree$tip.label
	
	# function to compute the states
	foo<-function(x){
	y<-sapply(x$maps,function(x) names(x)[1])
	names(y)<-x$edge[,1]
	y<-y[as.character(length(x$tip)+1:x$Nnode)]
	return(y)
	}
	
	
	#1. History of Symbiosis
	
	mydata<-read.csv("symbiont_charmatrix.txt", stringsAsFactors=FALSE,header=F,row.names=1)
	rownames(mydata)<-unlist(strsplit(rownames(mydata),split="_v1"))[c(T,F)]
	mydata
	
	rownames(mydata)
	mydata[tree$tip.label,]
	
	tree$states<-mydata$V2
	names(tree$states)<-rownames(mydata)
	symb<-tree$states
	mtrees<-make.simmap(tree,tree$states,model="ER", nsim=100)
	
	colsA<-colshcl[c(5,1)]; names(colsA)<-c(0,1)
	#plotSimmap(mtrees[[1]],cols,pts=F,ftype="off")
	
	AA<-sapply(mtrees,foo)
	piesA<-t(apply(AA,1,function(x,levels,Nsim) summary(factor(x,levels))/Nsim,levels=levels(factor(tree$states)), Nsim=100))
	plot.phylo(tree, use.edge.length = FALSE, no.margin=FALSE,show.tip.label=T, cex=1.3, main="Intracellular Eukaryotic Symbiont", node.depth = 0.4, align.tip.label = 2)
	nodelabels(pie=piesA,cex=1.2,piecol=colsA, cex=1.2)

	
	#2. History of Coloniality
	
	mydata2<-read.csv("colonial_char.txt", stringsAsFactors=FALSE,header=F,row.names=1)
	rownames(mydata2)<-unlist(strsplit(rownames(mydata2),split="_v1"))[c(T,F)]
	mydata2
	
	rownames(mydata2)
	mydata2[tree$tip.label,]
	
	tree$states<-mydata2$V2
	names(tree$states)<-rownames(mydata2)
		colonial<-tree$states

	mtrees<-make.simmap(tree,tree$states,model="ER", nsim=100)
	
	colsB<-colshcl[c(5,2)]; names(colsB)<-c(0,1)
	#plotSimmap(mtrees[[1]],cols,pts=F,ftype="off")
	
	AA<-sapply(mtrees,foo)
	piesB<-t(apply(AA,1,function(x,levels,Nsim) summary(factor(x,levels))/Nsim,levels=levels(factor(tree$states)), Nsim=100))
	plot.phylo(tree, use.edge.length = FALSE, no.margin=FALSE,show.tip.label=T, cex=1.3, main="Coloniality", node.depth = 0.4, adj = 0.1, align.tip.label = TRUE)
	nodelabels(pie=piesB,cex=1.2,piecol=colsB, cex=1.2)
	

	#3. History of Medusa
	
	mydata3<-read.csv("medusa_char.txt", stringsAsFactors=FALSE,header=F,row.names=1)
	rownames(mydata3)<-unlist(strsplit(rownames(mydata3),split="_v1"))[c(T,F)]
	mydata3
	
	rownames(mydata3)
	mydata3[tree$tip.label,]
	
	tree$states<-mydata3$V2
	names(tree$states)<-rownames(mydata3)
		medusa<-tree$states

	mtrees<-make.simmap(tree,tree$states,model="ER", nsim=100)
	
	colsC<-colshcl[c(5,3)]; names(colsC)<-c(0,1)
	#plotSimmap(mtrees[[1]],cols,pts=F,ftype="off")
	
	AA<-sapply(mtrees,foo)
	piesC<-t(apply(AA,1,function(x,levels,Nsim) summary(factor(x,levels))/Nsim,levels=levels(factor(tree$states)), Nsim=100))
	plot.phylo(tree, use.edge.length = FALSE, no.margin=FALSE,show.tip.label=T, cex=1.3, main="Medusa", node.depth = 0.4, adj = 0.1, align.tip.label = TRUE)
	nodelabels(pie=piesC,cex=1.2,piecol=colsC, cex=1.2)
	

	
	#4. History of Sessile Polyp
	
	mydata4<-read.csv("polyp_charmatrix.txt", stringsAsFactors=FALSE,header=F,row.names=1)
	rownames(mydata4)<-unlist(strsplit(rownames(mydata4),split="_v1"))[c(T,F)]
	mydata4
	
	rownames(mydata4)
	mydata4[tree$tip.label,]
	
	tree$states<-mydata4$V2
	names(tree$states)<-rownames(mydata4)
		polyp<-tree$states

	mtrees<-make.simmap(tree,tree$states,model="ER", nsim=100)
	
	colsD<-colshcl[c(5,4)]; names(colsD)<-c(0,1)
	#plotSimmap(mtrees[[1]],cols,pts=F,ftype="off")
		
	AA<-sapply(mtrees,foo)
	piesD<-t(apply(AA,1,function(x,levels,Nsim) summary(factor(x,levels))/Nsim,levels=levels(factor(tree$states)), Nsim=100))
	plot.phylo(tree, use.edge.length = FALSE, no.margin=FALSE,show.tip.label=T, cex=1.3, main="Polyp", node.depth = 0.4, adj = 0.1, align.tip.label = TRUE)
	nodelabels(pie=piesD,cex=1.2,piecol=colsD, cex=1.2)


	#plot all 4 characters and ASR on one tree	
plot.phylo(tree, no.margin=F,show.tip.label=T, use.edge.length=F, cex=1, main="Symbiosis, Coloniality, Medusa, Polyp", node.depth = 0.4, label.offset=2.5)
#nodelabels(text="          ",col="gray80",bg="gray80", cex=2.0,)
nodelabels(pch=16,col="gray90",bg="gray90", cex=4)

tips1<-which(symb[tree$tip.label]==1); tips0<-which(symb[tree$tip.label]==0)
tiplabels(tip=tips1,pch=16,col=colsA[2], cex=1.2, adj=1)
tiplabels(tip=tips0,pch=1,col=colsA[2], cex=1.2, adj=1)
nodelabels(pie=piesA,cex=.4,adj=c(.3,.9),piecol=colsA)
	
tips1<-which(colonial[tree$tip.label]==1); tips0<-which(colonial[tree$tip.label]==0)
tiplabels(tip=tips1,pch=16,col=colsB[2], cex=1.2, adj=1.5)
tiplabels(tip=tips0,pch=1,col=colsB[2], cex=1.2, adj=1.5)
nodelabels(pie=piesB,cex=.4,adj=c(.7,.9),piecol=colsB)
	
tips1<-which(medusa[tree$tip.label]==1); tips0<-which(medusa[tree$tip.label]==0)
tiplabels(tip=tips1,pch=16,col=colsC[2], cex=1.2, adj=2.0)
tiplabels(tip=tips0,pch=1,col=colsC[2], cex=1.2, adj=2.0)
nodelabels(pie=piesC,cex=.4,adj=c(.3,.1),piecol=colsC)

tips1<-which(polyp[tree$tip.label]==1); tips0<-which(polyp[tree$tip.label]==0)
tiplabels(tip=tips1,pch=16,col=colsD[2], cex=1.2, adj=2.5)
tiplabels(tip=tips0,pch=1,col=colsD[2], cex=1.2, adj=2.5)
nodelabels(pie=piesD,cex=.4,adj=c(.7,.1),piecol=colsD)
