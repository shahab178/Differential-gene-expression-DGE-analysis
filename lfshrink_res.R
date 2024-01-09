res_Syl_Vac <- lfcShrink(dds_Syl_Vac, coef = 2, type = "apeglm")
resOrdered <- res_Syl_Vac[order(res_Syl_Vac$padj),]
resSig <- subset(resOrdered, padj < 0.05)
lfshrink <- write.table(as.data.frame(resSig),"C:/Users/shahab/Desktop/tra/figado/lfshrink_DEGS.txt")


res <- results(dds_Syl_Vac)
resOrdered_res <- res[order(res$padj),]
resSig_res <- subset(resOrdered_res, padj < 0.05)
res <- write.table(as.data.frame(resSig_res),"C:/Users/shahab/Desktop/tra/figado/res_DEGS_res.txt")