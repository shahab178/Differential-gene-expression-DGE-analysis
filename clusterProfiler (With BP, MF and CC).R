#Import your dataset
#Col1 (Genes ID in Entrez form)
#Col2 (LogFoldChange)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("ReactomePA")
BiocManager::install("DOSE")
BiocManager::install("enrichplot")
BiocManager::install("ggnewscale")
BiocManager::install("reactome.db")
install.packages("stringr")
#Access libraries
library(DOSE)
library(enrichplot)
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggnewscale)
library(DOSE)
library(reactome.db)
library(rlang)
library(viridis)
library(ggplot2)
library(stringr)
#Run Enrichment

#geneList1 <- read.csv(/.../geneList1.txt", sep="")

d <- geneList
geneList <- d[,2]
valores <- d[,2]
names(geneList) <- as.character(d[,1])
genevector <- names(geneList)
geneList <- sort(geneList, decreasing = TRUE)

edo <- enrichPathway(genevector)
edo <- pairwise_termsim(edo, method = "JC", semData = NULL, showCategory = 200)
edo2 <- gsePathway(geneList, pvalueCutoff = 0.2)

egoCC <- enrichGO(genevector, ont="CC",  OrgDb = "org.Hs.eg.db")
egoCC <- pairwise_termsim(egoCC, method = "JC", semData = NULL, showCategory = 200)
egoCC2 <- gseGO(geneList, ont="CC",  OrgDb = "org.Hs.eg.db")

egoBP <- enrichGO(genevector, ont="BP",  OrgDb = "org.Hs.eg.db")
egoBP <- pairwise_termsim(egoBP, method = "JC", semData = NULL, showCategory = 200)
egoBP2 <- gseGO(geneList, ont="BP",  OrgDb = "org.Hs.eg.db")

egoMF <- enrichGO(genevector, ont="MF",  OrgDb = "org.Hs.eg.db")
egoMF <- pairwise_termsim(egoMF, method = "JC", semData = NULL, showCategory = 200)
egoMF2 <- gseGO(geneList, ont="MF",  OrgDb = "org.Hs.eg.db")

edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
egoxCC <- setReadable(egoCC, 'org.Hs.eg.db', 'ENTREZID')
egoxBP <- setReadable(egoBP, 'org.Hs.eg.db', 'ENTREZID')
egoxMF <- setReadable(egoMF, 'org.Hs.eg.db', 'ENTREZID')
##Temp
View(egoxBP@result)
goidBP <- c("GO:0006066", "GO:0006520", "GO:0006694", "GO:0031960")
egoxBP@result <- egoxBP@result[goidBP, ]
egoBP@result <- egoBP@result[goidBP, ]
cnetplot(egoxBP, foldChange=valores, circular = TRUE, colorEdge = TRUE, showCategory = 30)
dotplot(egoBP, showCategory=30) + ggtitle("dotplot for BP")
emapplot(egoBP, pie_scale=1.5,layout="kk")+scale_color_viridis()

p1 <- cnetplot(egoBP, foldChange=valores)
p2 <- cnetplot(egoBP, foldChange=valores, circular = TRUE, colorEdge = TRUE, showCategory = 5) + scale_color_viridis()
cowplot::plot_grid(p2, ncol=1, rel_widths=c(2, 2))
###############################################################
## PLOTS (ORA, Pahtways)

dotplot(edo, showCategory=20) + scale_color_viridis() + scale_y_discrete(labels=function(x) str_wrap(x, width=40)) + ggtitle("dotplot for ORA") 

#dotplot(edo, showCategory=20) + ggtitle("dotplot for ORA") + scale_color_viridis()
# Bar plot of enriched terms.
edo@result <- edo@result %>% arrange(desc(Count)) 

barplot(edo, showCategory=15, color="p.adjust")+ 
  ggtitle("barplot for ORA")+ 
  scale_color_viridis()
edo@result <- edo@result %>% arrange(desc(p.adjust))

#Ridgeline plot-Histogram
ridgeplot(edo2,showCategory=15, fill="pvalue")+ 
  #fill= "p.adjust" if necessary
  xlab("logFC")

#Extra
upsetplot(edo)+ ggtitle("upsetplot for ORA") +scale_color_viridis()
#Network plot of enriched terms.
p1 <- cnetplot(edox, foldChange=valores)
p2 <- cnetplot(edox, foldChange=valores, circular = TRUE, colorEdge = TRUE, showCategory = 5) + scale_color_viridis()
cowplot::plot_grid(p2, ncol=1, rel_widths=c(2, 2))

#Labelling nodes by selected subset
p3 <- cnetplot(edox, node_label="none")
p4 <- cnetplot(edox, node_label="gene") 
p5 <- cnetplot(edox, node_label="category") 
p6 <- cnetplot(edox, node_label="all") 
cowplot::plot_grid(p6, ncol=1)

#Heatmap plot of enriched terms
p7 <- heatplot(edo)
p8 <- heatplot(edo, foldChange=valores)+ scale_color_viridis()
cowplot::plot_grid(p8, ncol=1)

#Plot for results obtained from hypergeometric test and gene set enrichment analysis.
p9 <- emapplot(edo)
p10 <- emapplot(edo, cex_category=1.5)
p11 <- emapplot(edo,layout="circle")
p12 <- emapplot(edo, cex_category=1.5,layout="circle", showCategory = 20, cex_line = 0.4)+ scale_color_viridis()
cowplot::plot_grid(p12, ncol=1)

#Network using specific pathways
View(edox@result)
#choose your pathways
goid <- c("R-HSA-909733", "R-HSA-913531", "R-HSA-877300", "R-HSA-983170", "R-HSA-1236975", "R-HSA-198933")
edox@result <- edox@result[goid, ]
cnetplot(edox, foldChange=valores, circular = TRUE, colorEdge = TRUE)+ scale_color_viridis()
emapplot(edox, cex_category=1.5,layout="kk")+ scale_color_viridis() 
cnetplot(edox, foldChange=valores, colorEdge = TRUE)+ scale_color_viridis() 

#After running this, you should run the enrichment again
edo <- enrichPathway(genevector)
edo <- pairwise_termsim(edo, method = "JC", semData = NULL, showCategory = 200)
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
###########################################################
## PLOTS (CC)##
#Dot plot of enriched terms.
rm (list = c("p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "p10", "p11", "p12"))
dotplot(egoCC, showCategory=20) + scale_color_viridis() + scale_y_discrete(labels=function(x) str_wrap(x, width=40)) + ggtitle("dotplot for CC")
#dotplot(egoCC, showCategory=15) + ggtitle("dotplot for CC") + scale_color_viridis()
# Bar plot of enriched terms.
egoCC@result <- egoCC@result %>% arrange(desc(Count))

barplot(egoCC, showCategory=15, color="p.adjust")+ 
  ggtitle("barplot for CC")+ 
  scale_color_viridis()

egoCC@result <- egoCC@result %>% arrange(desc(p.adjust))
#Ridgeline plot-histogram
ridgeplot(egoCC2, fill="pvalue")+ 
  #fill= "p.adjust" if necessary
  xlab("logFC")
#Extra
upsetplot(egoCC)+ ggtitle("upsetplot for CC") +scale_color_viridis()
#Network plot of enriched terms.
p1 <- cnetplot(egoxCC, foldChange=valores)
p2 <- cnetplot(egoxCC, foldChange=valores, circular = TRUE, colorEdge = TRUE, showCategory = 5) + scale_color_viridis()
cowplot::plot_grid(p2, ncol=1, rel_widths=c(1, 1))

#Labelling nodes by selected subset
p3 <- cnetplot(egoxCC, node_label="none")
p4 <- cnetplot(egoxCC, node_label="gene") 
p5 <- cnetplot(egoxCC, node_label="category") 
p6 <- cnetplot(egoxCC, node_label="all") 
cowplot::plot_grid(p6, ncol=1)

#Heatmap plot of enriched terms
p7 <- heatplot(egoxCC)
p8 <- heatplot(egoxCC, foldChange=valores) + scale_color_viridis()
cowplot::plot_grid(p8, ncol=1)

#Plot for results obtained from hypergeometric test and gene set enrichment analysis.
p9 <- emapplot(egoCC)
p10 <- emapplot(egoCC, cex_category=1.5)
p11 <- emapplot(egoCC,layout="kk")
p12 <- emapplot(egoCC, cex_category=1.5,layout="kk")  + scale_color_viridis()
cowplot::plot_grid(p12, ncol=1)
#Network using specific pathways
View(egoxCC@result)
#choose your pathways
goidCC <- c("GO:1990204", "GO:0098563", "GO:0030665", "GO:0044815", "GO:0031985")
egoxCC@result <- egoxCC@result[goidCC, ]
cnetplot(egoxCC, foldChange=valores, circular = TRUE, colorEdge = TRUE)+ scale_color_viridis()
emapplot(egoxCC, cex_category=1.5,layout="kk")+ scale_color_viridis() 
cnetplot(egoxCC, foldChange=valores, colorEdge = TRUE)+ scale_color_viridis() 

#After running this, you should run the enrichment again
egoCC <- enrichGO(genevector, ont="CC",  OrgDb = "org.Hs.eg.db")
egoCC <- pairwise_termsim(egoCC, method = "JC", semData = NULL, showCategory = 200)
egoxCC <- setReadable(egoCC, 'org.Hs.eg.db', 'ENTREZID')


###############################################################
##PLOTS (BP)##
rm (list = c("p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "p10", "p11", "p12"))
#Dot plot of enriched terms.
dotplot(egoBP, showCategory=20) + scale_color_viridis() + scale_y_discrete(labels=function(x) str_wrap(x, width=40)) + ggtitle("dotplot for BP")
#dotplot(egoBP, showCategory=15) + ggtitle("dotplot for BP") + scale_color_viridis()

# Bar plot of enriched terms.
egoBP@result <- egoBP@result %>% arrange(desc(Count))

barplot(egoBP, showCategory=15, color="p.adjust")+ 
  ggtitle("barplot for BP")+ 
  scale_color_viridis()

egoBP@result <- egoBP@result %>% arrange(desc(p.adjust))

#Ridgeline plot-histogram
ridgeplot(egoBP2, fill="pvalue")+
  #fill= "p.adjust" if necessary
  xlab("logFC")
#Extra
upsetplot(egoBP)+ ggtitle("upsetplot for BP") +scale_color_viridis()

#Network plot of enriched terms.
#Layout of the map, e.g. 'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl' or 'lgl'.
p1 <- cnetplot(egoxBP, foldChange=valores)
p2 <- cnetplot(egoxBP, foldChange=valores, circular = TRUE, colorEdge = TRUE, showCategory = 5, layout = "star") + scale_color_viridis()
cowplot::plot_grid(p2, ncol=1, rel_widths=c(1, 1)) 

#Labelling nodes by selected subset
p3 <- cnetplot(egoxBP, node_label="none")
p4 <- cnetplot(egoxBP, node_label="gene") 
p5 <- cnetplot(egoxBP, node_label="category") 
p6 <- cnetplot(egoxBP, node_label="all") 
cowplot::plot_grid(p6, ncol=1)

#Heatmap plot of enriched terms
p7 <- heatplot(egoxBP)
p8 <- heatplot(egoxBP, foldChange=valores) + scale_color_viridis()
cowplot::plot_grid(p8, ncol=1)

#Plot for results obtained from hypergeometric test and gene set enrichment analysis.
p9 <- emapplot(egoBP)
p10 <- emapplot(egoBP, cex_category=1.5)
p11 <- emapplot(egoBP,layout="kk")
p12 <- emapplot(egoBP, cex_category=1.5,layout="kk")+ scale_color_viridis() 
cowplot::plot_grid(p12, ncol=1)

#Network using specific pathways.Change according the refs below
View(egoxBP@result)

#choose your pathways. change according the refs below
goidBP <- c("GO:0060337", "GO:0071357", "GO:0034340", "GO:0060333")
egoxBP@result <- egoxBP@result[goidBP, ]
cnetplot(egoxBP, foldChange=valores, circular = TRUE, colorEdge = TRUE)+ scale_color_viridis()
emapplot(egoxBP, cex_category=1.5,layout="kk")+ scale_color_viridis() 
cnetplot(egoxBP, foldChange=valores, colorEdge = TRUE)+ scale_color_viridis()

barplot(egoxBP, showCategory=15, color="p.adjust")+ 
  ggtitle("barplot for BP")+ 
  scale_color_viridis()
dotplot(egoxBP, showCategory=20) + ggtitle("dotplot for BP")+ scale_color_viridis()
#Ridgeline plot - Histogram
ridgeplot(egoBP2, fill="pvalue")+
  #fill= "p.adjust" if necessary
  xlab("logFC")
heatplot(egoxBP, foldChange=valores)+ scale_color_viridis()

#After running this, you should run the enrichment again
egoBP <- enrichGO(genevector, ont="BP",  OrgDb = "org.Hs.eg.db")
egoBP <- pairwise_termsim(egoBP, method = "JC", semData = NULL, showCategory = 200)
egoxBP <- setReadable(egoBP, 'org.Hs.eg.db', 'ENTREZID')
###############################################################
##PLOTS (MF)##
rm (list = c("p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "p10", "p11", "p12"))
#Dot plot of enriched terms.
dotplot(egoMF, showCategory=20) + scale_color_viridis() + scale_y_discrete(labels=function(x) str_wrap(x, width=40)) + ggtitle("dotplot for MF")
#dotplot(egoMF, showCategory=15) + ggtitle("dotplot for MF") + scale_color_viridis()

# Bar plot of enriched terms.
egoMF@result <- egoMF@result %>% arrange(desc(Count))

barplot(egoMF, showCategory=15, color="p.adjust")+ 
  ggtitle("barplot for MF")+ 
  scale_color_viridis()

egoMF@result <- egoMF@result %>% arrange(desc(p.adjust))

#Ridgeline plot-Histogram
ridgeplot(egoMF2, showCategory = 20, fill="pvalue") + scale_y_discrete(labels=function(x) str_wrap(x, width=40)) 
  #fill= "p.adjust" if necessary
  xlab("logFC")
#Extra
upsetplot(egoMF)+ ggtitle("upsetplot for MF") +scale_color_viridis()
#Network plot of enriched terms.
p1 <- cnetplot(egoxMF, foldChange=valores)
p2 <- cnetplot(egoxMF, foldChange=valores, circular = TRUE, colorEdge = TRUE, showCategory = 5) + scale_color_viridis()
cowplot::plot_grid(p2, ncol=1, rel_widths=c(1, 1))
#Labelling nodes by selected subset
p3 <- cnetplot(egoxMF, node_label="none")
p4 <- cnetplot(egoxMF, node_label="gene") 
p5 <- cnetplot(egoxMF, node_label="category") 
p6 <- cnetplot(egoxMF, node_label="all") 
cowplot::plot_grid(p6, ncol=1)

#Heatmap plot of enriched terms
p7 <- heatplot(egoxMF)
p8 <- heatplot(egoxMF, foldChange=valores) + scale_color_viridis()
cowplot::plot_grid(p8, ncol=1)

#Plot for results obtained from hypergeometric test and gene set enrichment analysis.
p9 <- emapplot(egoMF)
p10 <- emapplot(egoMF, cex_category=1.5)
p11 <- emapplot(egoMF,layout="kk")
p12 <- emapplot(egoMF, cex_category=1.5,layout="kk")+ scale_color_viridis() 
cowplot::plot_grid(p12, ncol=1)

#Network using specific pathways. change according the refs below
View(egoxMF@result)
#choose your pathways
goidMF <- c("GO:0016887", "GO:0003678", "GO:0003777", "GO:0003688", "GO:0003774")
egoxMF@result <- egoxMF@result[goidMF, ]
cnetplot(egoxMF, foldChange=valores, circular = TRUE, colorEdge = TRUE)+ scale_color_viridis()
emapplot(egoxMF, cex_category=1.5,layout="kk")+ scale_color_viridis() 
cnetplot(egoxMF, foldChange=valores, colorEdge = TRUE)+ scale_color_viridis() 

#After running this, you should run the enrichment again
egoMF <- enrichGO(genevector, ont="MF",  OrgDb = "org.Hs.eg.db")
egoMF <- pairwise_termsim(egoMF, method = "JC", semData = NULL, showCategory = 200)

egoxMF <- setReadable(egoMF, 'org.Hs.eg.db', 'ENTREZID')
#######################################################################
#http://yulab-smu.top/clusterProfiler-book/chapter4.html#gsencg-fuction
#Density plots for diseases association with the list of DEGs

#gseDO means GSE Diseases Ontology
y <- gseDO(geneList,
           nPerm         = 100,
           minGSSize     = 120,
           pvalueCutoff  = 1,
           verbose       = FALSE)
ridgeplot(y, fill= "pvalue")

#gseDGN means GSE Diseases Ontology

dgn <- gseDGN(geneList,
              nPerm         = 100,
              minGSSize     = 120,
              pvalueCutoff  = 1,
              verbose       = FALSE)

ridgeplot(dgn, fill= "pvalue")

#gseNCG means GS Network of Cancer Gene

ncg <- gseNCG(geneList,
              nPerm         = 100,
              minGSSize     = 120,
              pvalueCutoff  = 1,
              verbose       = FALSE)

ridgeplot(ncg, fill= "pvalue")
###tree plots
#https://www.bioconductor.org/packages/devel/bioc/manuals/enrichplot/man/enrichplot.pdf
#https://rdrr.io/github/GuangchuangYu/enrichplot/man/treeplot.html
rm (list = c("p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "p10", "p11", "p12"))
#ORA#
p1 <- treeplot(edox, offset_tiplab = rel(0.5)) + scale_color_viridis() + theme(legend.position = "bottom")
p2 <- treeplot(edox, hclust_method = "average", offset_tiplab = rel(0.5)) + scale_color_viridis() + theme(legend.position = "bottom")
aplot::plot_list(p1, p2, tag_levels='A')

p2 <- treeplot(edox, hclust_method = "average", offset_tiplab = rel(0.5), showCategory = 50) + scale_color_viridis() + theme(legend.position = "bottom")
aplot::plot_list(p2)
#CC#
p4 <- treeplot(egoxCC, hclust_method = "average", offset_tiplab = rel(0.5), showCategory = 50) + scale_color_viridis() + theme(legend.position = "bottom")
aplot::plot_list(p4)
#MF#
p6 <- treeplot(egoxMF, hclust_method = "average", offset_tiplab = rel(0.5), showCategory = 50) + scale_color_viridis() + theme(legend.position = "bottom")
aplot::plot_list(p6)
#BP#
p8 <- treeplot(egoxBP, hclust_method = "average", offset_tiplab = rel(0.5), showCategory = 50) + scale_color_viridis() + theme(legend.position = "bottom")
aplot::plot_list(p8)
##################################################################
#Export tables

install.packages("writexl")
library(writexl)

tabelaORA <- edox@result
tabelaCC <- egoxCC@result
tabelaBP <- egoxBP@result
tabelaMF <- egoxMF@result


write_xlsx(tabelaORA, "/.../tabelaORA.xlsx")
write_xlsx(tabelaCC, "/.../tabelaCC.xlsx")
write_xlsx(tabelaBP, "/.../tabelaBP.xlsx")
write_xlsx(tabelaMF, "/.../tabelaMF.xlsx")
