library(ggplot2)

data <- read.table('mat.txt')
data
#plot occ vs tax, size = # parts
p <- ggplot(data, aes(occ, tax))
p + geom_point()
p + geom_point(aes(size = parts)) 

#plot occ vs parts, size = # parts
p <- ggplot(data, aes(occ, parts))
p + geom_point()
p + geom_point(aes(size = tax)) 