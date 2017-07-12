
library(ape)
#library(devtools)
#install_github("vnminin/indorigin",build_vignettes=TRUE)

library(indorigin)

#install.packages("~/Downloads/indorigin_0.0.1.tar.gz",type="source",repos=NULL)
#install.packages("indorigin", repos="http://R-Forge.R-project.org")

trees = read.tree("PhyloBayesCATGTRI_OF75tx_2.treelist.pruned")
#trees = read.tree("testPPtrees")
Trees<-root(trees,outgroup="Mnemiopsis_leidyi_v1_Mnemiopsis_leidyi",resolve.root=T)
length(Trees[[1]]$tip.label)



runIndOrigins<-function(myRootedTrees,myTraits,numbers2test,prGain,prLoss)
{	save<-c()
	for (i in 1:length(num_origins))
	{
	sabrinaIndOriginResultsM4 = testIndOrigin(inputTrees=myRootedTrees, traitData=myTraits, initLambda01=.01, initLambda10=0.1, priorAlpha01=1,
 		priorBeta01=prGain, priorAlpha10=1, priorBeta10=prLoss, mcmcSize=55000, mcmcBurnin=5000, mcmcSubsample=5, mcSize=50000, testThreshold=numbers2test[i])
	print(paste("BF for number origins =<",num_origins[i]))
	BFs<-getBF(sabrinaIndOriginResultsM4)
	prior<-getPriorProb(sabrinaIndOriginResultsM4)[1]
	post<-getPostProb(sabrinaIndOriginResultsM4)[1]
	save<-c(save,numbers2test[i],prior, post, BFs)
	}
	return(as.matrix(save,ncol=length(numbers2test)))
}


print("coloniality")
Traits = read.csv("colonial_charmatrix.txt", header=FALSE)

tiptrait.num = dim(Traits)[1]; TraitVec = numeric(tiptrait.num)
tipNames = as.character(Traits[,1]); TraitVec = Traits[,2]; names(TraitVec) = tipNames
num_origins<-c(0,1,2,3,4)
run1<-runIndOrigins(Trees,TraitVec,num_origins,1,10)
mat<-matrix(c(run1),ncol=5)
print(mat)

print("polyps")
Traits = read.csv("polyp_charmatrix.txt", header=FALSE)

tiptrait.num = dim(Traits)[1]; TraitVec = numeric(tiptrait.num)
tipNames = as.character(Traits[,1]); TraitVec = Traits[,2]; names(TraitVec) = tipNames
num_origins<-c(0,1,2,3,4)
run1<-runIndOrigins(Trees,TraitVec,num_origins,1,10)
mat<-matrix(c(run1),ncol=5)
print(mat)

print("medusa")
Traits = read.csv("medusa_char.txt", header=FALSE)

tiptrait.num = dim(Traits)[1]; TraitVec = numeric(tiptrait.num)
tipNames = as.character(Traits[,1]); TraitVec = Traits[,2]; names(TraitVec) = tipNames
num_origins<-c(0,1,2,3,4)
run1<-runIndOrigins(Trees,TraitVec,num_origins,1,10)
mat<-matrix(c(run1),ncol=5)
print(mat)


print("symbionts")
Traits = read.csv("symbiont_charmatrix.txt", header=FALSE)

tiptrait.num = dim(Traits)[1]; TraitVec = numeric(tiptrait.num)
tipNames = as.character(Traits[,1]); TraitVec = Traits[,2]; names(TraitVec) = tipNames
num_origins<-c(0,1,2,3,4)
run1<-runIndOrigins(Trees,TraitVec,num_origins,1,10)
mat<-matrix(c(run1),ncol=5)
print(mat)
