if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')


library(EnhancedVolcano)
library(ggplot2)

##################################################################
# USING P.ADJ
#EnhancedVolcano <- read.delim("/.../EnhancedVolcano.txt", row.names=1)

keyvals <- ifelse(
  EnhancedVolcano$log2FoldChange < 1 & EnhancedVolcano$padj < 0.05, 'royalblue',
  ifelse(EnhancedVolcano$log2FoldChange > 1 & EnhancedVolcano$padj < 0.05, 'red',
         'black'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'Up-regulated'
names(keyvals)[keyvals == 'black'] <- 'NonDEGs'
names(keyvals)[keyvals == 'royalblue'] <- 'Down-regulated'

EnhancedVolcano(EnhancedVolcano,
                lab =rownames(EnhancedVolcano),
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-10,10),
                ylim = c(-0.05, -log10(10e-8)),
                xlab = bquote(~Log[2] ~ "fold change"),
                ylab = bquote(~-Log[10] ~ italic(p)~"adj"),
                pointSize = 2,
                selectLab = rownames(EnhancedVolcano)[which(names(keyvals) %in% c('Up-regulated', 'Down-regulated'))],
                labSize = 3,
                cutoffLineType = 'twodash',
                title = "Liver",
                subtitle = "Sylvativ Vs Control",
                legendPosition = "right",
                legendLabSize = 9,
                colAlpha = 0.7,
                labCol = 'black',
                border = 'full',
                borderColour = 'black',
                borderWidth = 1,
                drawConnectors = FALSE,
                cutoffLineCol = 'black',
                cutoffLineWidth = 0.8,
                colCustom= keyvals,
                gridlines.major = FALSE,
                gridlines.minor = FALSE
) +
  ggplot2::annotate("label", x =c(-9.3,0,9.3), y = 6.0 , label = c("2171","29403","1873"), col=c("royalblue", "black", "red"))

# change the quantity of Down-regulated, Up-regulated and  the total genes manually in the last command


