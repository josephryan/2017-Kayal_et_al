setwd(source("~/Dropbox/Dave_shared/map gene partitions on tree/partition overlap comparison.R")

myx<-c("Kudoa_iwatai" ,"Myxobolus_cerebralis","Myxobolus_pendula"  ,"Thelohanellus_kitauei")

noncnid<-c("Trichoplax_adhaerens","Strongylocentrotus_purpuratus" ,"Taeniopygia_guttata" ,"Amphimedon_queenslandica","Drosophila_melanogaster","Capitella_teleta","Lottia_gigantea","Mnemiopsis_leidyi")

par(mfrow=c(2,2),oma=c(1,1,5,1))


#Quantify partition overlap for OF-PTP
datof<-read.table("OF50partitiontable.txt",row.names=1,stringsAsFactors=F)
rows<-rownames(datof)
datof[datof =="white"]<-0; datof[datof =="black"]<-1
datof <-data.frame(lapply(datof,as.numeric));rownames(datof)<-rows

cnid<-rownames(datof)[-(which(rownames(datof) %in% c(myx,noncnid)))]
myxSum1<-names(which(colSums(datof[myx,])>=1))  #OG groups containing Myxo taxa
myxSum2<-names(which(colSums(datof[myx,])>=2))  #OG groups containing Myxo taxa
myxSum3<-names(which(colSums(datof[myx,])>=3))  #OG groups containing Myxo taxa
myxSum4<-names(which(colSums(datof[myx,])>=4))  #OG groups containing Myxo taxa

write(myxSum1,file="OF_m1")
write(myxSum2,file="OF_m2")
write(myxSum3,file="OF_m3")
write(myxSum4,file="OF_m4")

cnidSum <-names(which(colSums(datof[cnid,])>1))
otherSum<-names(which(colSums(datof[noncnid,])>1))

myxProp<-colSums(datof[myx,])/length(myx)
cnidProp<-colSums(datof[cnid,])/length(cnid)
otherProp<-colSums(datof[noncnid,])/length(noncnid)
venn.input<-list(Myxozoa=myxSum, Cnidaria=cnidSum, others=otherSum)
venn(venn.input)
title("OrthoFinder data partitions shared by >=1 member",outer=T)
hist(myxProp[which(colSums(datof[myx,])>=1)], main="Myx")
hist(cnidProp[which(colSums(datof[cnid,])>=1)],main="cnidaria")
hist(otherProp[which(colSums(datof[noncnid,])>=1)],main="outgroups")


#Quantify partition overlap for Agalma

datag<-read.table("agalma.partitiontable.txt",row.names=1,stringsAsFactors=F)
rows<-rownames(datag)
datag[datag=="white"]<-0; datag[datag=="black"]<-1
datag <-data.frame(lapply(datag,as.numeric)); rownames(datag)<-rows

cnid<-rownames(datag)[-(which(rownames(datag) %in% c(myx,noncnid)))]
myxSum<-names(which(colSums(datag[myx,])>=1))  #Agalma groups containing Myxo taxa
myxSum1<-names(which(colSums(datof[myx,])>=1))  #OG groups containing Myxo taxa
myxSum2<-names(which(colSums(datof[myx,])>=2))  #OG groups containing Myxo taxa
myxSum3<-names(which(colSums(datof[myx,])>=3))  #OG groups containing Myxo taxa
myxSum4<-names(which(colSums(datof[myx,])>=4))  #OG groups containing Myxo taxa

write(myxSum1,file="AG_m1")
write(myxSum2,file="AG_m2")
write(myxSum3,file="AG_m3")
write(myxSum4,file="AG_m4")



cnidSum <-names(which(colSums(datag[cnid,])>1))
otherSum<-names(which(colSums(datag[noncnid,])>1))

myxProp<-colSums(datag[myx,])/length(myx)
cnidProp<-colSums(datag[cnid,])/length(cnid)
otherProp<-colSums(datag[noncnid,])/length(noncnid)



venn.input<-list(Myxozoa=myxSum, Cnidaria=cnidSum, others=otherSum)
venn(venn.input)
title("Agalma data partitions shared by >=1 member",outer=T)

hist(myxProp[which(colSums(datag[myx,])>=1)], main="Myx")
hist(cnidProp[which(colSums(datag[cnid,])>=1)],main="cnidaria")
hist(otherProp[which(colSums(datag[noncnid,])>=1)],main="outgroups")



