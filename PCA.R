#FROM: 
#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/
#https://f0nzie.github.io/machine_learning_compilation/detailed-study-of-principal-component-analysis.html
#https://rstudio-pubs-static.s3.amazonaws.com/323416_ab58ad22d9e64ba2831569cf3d14a609.html
#http://www.sthda.com/english/wiki/fviz-pca-quick-principal-component-analysis-data-visualization-r-software-and-data-mining

install.packages("factoextra")
install.packages("ggExtra")  #from https://www.rdocumentation.org/packages/ggExtra/versions/0.9/topics/ggMarginal


library("factoextra")
library("ggExtra")
library("ggplot2")
library("scales")

groups <- as.factor(PCAf$Group)

PCAf[1] <- NULL  

res.pca <- prcomp(PCAf, scale = TRUE)

#Biplot sem os vetores
p <- fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             geom = c("point"),
             palette = c("#eaac43", "#ff7373", "#1657e1", "#38c3e1"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE,
             pointsize = 3
) + theme(legend.position = "bottom") + ggtitle("PCA") 
print(p)

#Adding histogram
p1 <- ggExtra::ggMarginal(p, type = "histogram", groupFill = TRUE, groupColour = TRUE)
print(p1)

#BIPLOT com os vetores dos genes
p <- fviz_pca_biplot(res.pca,
                  col.ind = groups, # color by groups,
                  geom = c("point"),
                  col.var = "grey",
                  alpha.var = 0.2,
                  palette = c("#eaac43", "#ff7373", "#1657e1", "#38c3e1"),
                  addEllipses = TRUE, # Concentration ellipses
                  ellipse.type = "confidence",
                  legend.title = "Groups: ",
                  repel = TRUE,
                  pointsize=3
) + theme(legend.position = "bottom") + ggtitle("PCA") 
print(p)

#Adding histogram
p1 <- ggExtra::ggMarginal(p, type = "histogram", groupFill = TRUE, groupColour = TRUE)
print(p1)



#to change the theme see: https://ggplot2.tidyverse.org/reference/theme.html

#if need to ajust the x and y axis limits add after the parenteses: + xlim(-1, 15) + ylim(-5, 5)

#For additional color see https://www.color-hex.com/

