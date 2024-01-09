# Differential gene expression DGE analysis

## OBS:

1. /.../ should be change to your convinient path in each system.
2. Pre processing and Data Quality pipeline can be found at [NGS](https://github.com/shahab178/NGS_illumina).

## Differential Expression Analysis:

For differential expression analysis, several software packages and pipelines as well as
numerous Bioinformatics programs for RNA Seq data analysis using different algorithms and
models have been developed. The diversity of the methodology makes it possible to customize
analysis protocols (i.e., pipelines) by choosing and concatenating programs that provide the
best fit for each specific goal and situation. We have been using DESeq2 for the count based 
statistics assuming a negative binomial distribution model.
[kallisto_tximport_Deseq2.R](https://github.com/shahab178/Differential-gene-expression-analysis/blob/main/kallisto_tximport_Deseq2.R)

## Exploratory Analysis and Quality Control
- PCA
- Volcano plot with EnhancedVolcano package
- ComplexHeatmap
- ClusterProfiler package for Gene Ontology (GO)
- ifshrinkage test (for Log2FC normalization)
Generally, visualization of RNA seq data is similar to any other type of genomic
sequencing data, which can be done at the level of reads, total counts, normalized or unnormalized,
cluster visualization, principal component analysis (PCA), and heatmaps.
Exploratory analysis is essential for data understanding and quality control.
It can provide us a sense of the most obvious patterns in the data and assist
us in identifying quality issues, sample swaps, and contamination. For this
reduction of the dimension of the data set, PCA used the linear combination of the
original data (gene expression values) to define a new set of unrelated variables (principal components).
RNA sequencing with high throughput is a potent method for analyzing gene expression.
Extreme deviation of a sample from samples of the same treatment group may occur as a result
of technical variance or real biological variations due to the complex multi step methods in data collecting.
