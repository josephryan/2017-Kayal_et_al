# install.packages('VennDiagram')
library(VennDiagram)

draw.pairwise.venn(175, 376, 94, category = c("OPP 682", "Agamla 2812"), lty = rep("blank", 
    2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 
    0), cat.dist = rep(0.025, 2))