library(ggfortify)

#read table (RNA_seq)
df = read.csv("data/countMatrix_sorted_polyA.csv", row.name=1, sep=",", header = T, check.names = FALSE)
df$gene_name = NULL
df$gene_type = NULL

datTraits = read.csv("data/experimentList_sorted.csv", row.name=1, sep=",", header = T)

df_t <- t(df)
df_pca <- prcomp(df_t)


pdf("PCA.pdf", width = 9, height = 7)
p = autoplot(df_pca, data=datTraits, colour="Group", shape="Sex") + 
  theme_bw()
print(p)
dev.off()

