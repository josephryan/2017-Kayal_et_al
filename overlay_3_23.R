
#install.packages("reshape")
#install.packages("ggplot2")

library(ggplot2)
library(reshape2)
#setwd("~/Downloads")

data <- read.table("zap_sum", header = TRUE)
dim(data)
method<-rep("Zapata (1262 parts)",nrow(data))
data<-cbind(method,data)

PI_A <-data1$Parsimony_informative_sites
numb_taxa_B <-data$No_of_taxa
align_len_B <-data$Alignment_length
no_v_sites_B <-data$No_variable_sites

data1 <- read.table("AG_50_summary.txt", header = TRUE)
method<-rep("AGALMA_50 (962) parts)",nrow(data1))
dim(data1)
data1<-cbind(method,data1)

PI_B <-data1$Parsimony_informative_sites
numb_taxa_C <-data1$No_of_taxa
align_len_C <-data1$Alignment_length
no_v_sites_C <-data1$No_variable_sites


data2 <- read.table("chang_sum", header = TRUE)
dim(data2)
method<-rep("Chang (200 parts)",nrow(data2))
data2<-cbind(method,data2)

PI_C <-data2$Parsimony_informative_sites
numb_taxa_A <-data2$No_of_taxa
align_len_A <-data2$Alignment_length
no_v_sites_A <-data2$No_variable_sites


data3 <- read.table("OF_50_summary.txt", header = TRUE)
method<-rep("OF_50 (373 parts)",nrow(data3))
dim(data3)
data3<-cbind(method,data3)

PI_D <-data3$Parsimony_informative_sites
numb_taxa_D <-data3$No_of_taxa
align_len_D <-data3$Alignment_length
no_v_sites_D <-data3$No_variable_sites


#manual data reshape:
dat<-rbind(data,data1,data2,data3)
dat<-dat[,which(names(dat) %in% c("method","Parsimony_informative_sites","No_of_taxa","Alignment_length","No_variable_sites"))]

		##test
# hist(PI_A)
# da<-density(PI_A);db<-density(PI_B);dc<-density(PI_C)
# plot(da,xlab="PI sites?",ylab="Density");lines(db,col="red");lines(dc,col="blue")
# legend("topright",legend=c("OF50","Borowiec","CEGMA"),text.col=c("black","red","blue"))


## grid graphics (as opposed to traditional/default graphics) are prone to glitches
## if below fails, try: Window>New Quartz Device Window (shift-cmd-N)

p1<-qplot(Parsimony_informative_sites, data=dat, geom="histogram", fill= method, alpha=I(0.5), main="Distribution of Parsimony Informative Sites", xlab="PI sites", ylab="Density")

p2<-qplot(No_of_taxa, data=dat, geom="histogram", fill= method, alpha=I(0.5), main="Distribution of Number of Taxa", xlab="Number of Taxa", ylab="Density")
      
p3<-qplot(Alignment_length, data=dat, geom="histogram", fill= method, alpha=I(0.5), main="Distribution of Alignment Lengths", xlab="Alignment Lengths", ylab="Number of Sequences")
      
p4<-qplot(No_variable_sites, data=dat, geom="histogram", fill= method, alpha=I(0.5), main="Distribution of Number of Variable Sites", xlab="Number of Variable Sites", ylab="Density")

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#pdf("densities2.pdf",width=11,height=8)
multiplot(p1,p2,p3,p4,cols=2)
#dev.off()