
#install.packages("circlize")
#BiocManager::install("ComplexHeatmap")
library("ComplexHeatmap")
library("circlize")
library("writexl")

#convert the dataframe in matrix
heatEX <-as.matrix(heatEX)

ha1 = HeatmapAnnotation(
  dist1 = anno_barplot(
    colSums(heatEX), 
    bar_width = 1, 
    gp = gpar(col = "blue", fill = "orange"), 
    border = FALSE),
  show_annotation_name = FALSE)
ha2 = rowAnnotation(
  dist2 = anno_barplot(
    rowSums(heatEX), 
    bar_width = 1, 
    gp = gpar(col = "blue", fill = "orange"), 
    border = FALSE
  ), show_annotation_name = FALSE)

col_fun = RColorBrewer::brewer.pal(n =9, name = "GnBu") 
#to see more colors: https://www.r-graph-gallery.com/38-rcolorbrewers-palettes.html

ht_list = Heatmap(heatEX, name = "expression",
                  cluster_columns = TRUE, show_row_dend = TRUE, #rect_gp = gpar(col= "white"), 
                  show_column_names = TRUE,  show_row_names = FALSE, column_names_gp = gpar(fontsize = 6),
                  row_names_side = "left", row_names_gp = gpar(fontsize = 9),
                  #column_names_side = "top",
                  column_title = 'Gene Expression',
                  top_annotation = ha1, col = col_fun) + ha2
                   
draw(ht_list, ht_gap = unit(1, "mm"))

