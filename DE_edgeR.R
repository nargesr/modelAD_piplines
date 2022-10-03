setwd("I:/Lab/Model AD/Differential Gene Expression/")
library(edgeR)
library(ggplot2)
library("ggrepel")

###---------------------- 
# Read data
#read table (RNA_seq or ATAC_seq)
df = read.csv("Data/expressionList_3xTgAD_hipp_female_18mon", sep=",", header = T)

#read gene list to convert gene ID to gene list
annot = read.table("Data/geneList", sep="\t", header = T)

row.names(df) = df$gene_id
df$gene_id = NULL

## convert gene ID to gene Name
for (i in c(53496:dim(df)[1])){
  row.names(df)[i] = strsplit(row.names(df)[i], "\\.")[[1]][1]
  if (length(which(row.names(df)[i] == annot$Gene.ID)) == 1) {
    row.names(df)[i] = annot$Gene.Symbol[which(row.names(df)[i] == annot$Gene.ID)]
  }
}

# create group by header
groups = rep(0, dim(df)[2])
for (i in c(1:dim(df)[2])) {
  tmp = strsplit(colnames(df)[i], "_")[[1]]
  if (tmp[3] == "3xTgADWT"){
    groups[i] = 1
  }
  if (tmp[3] == "3xTgADHO"){
    groups[i] = 2
  }
}

###---------------------- 
## ANALYSING
#change the group numbers (group,#samples)
#We always want our control to be group1
#because when we see the results is always if it's up or down in group2
#For example if my controls are group 2 and my 5xFAD are group1
#When I see the results from the et table, such as logFC, will be related to 5xFAD (group2)

y = DGEList(counts = df, group=groups)
y = calcNormFactors(y)
design = model.matrix(~groups)
y = estimateDisp(y, design)

keep = filterByExpr(y)
y = y[keep, , keep.lib.sizes=FALSE]

#Testing for DE Genes
et = exactTest(y)

#extract table from the exact test( here is where we know if they are DE or not)
et_out = (topTags(et, n=Inf, adjust.method = "BH"))
et = et_out$table
et$gene_name = rownames(et)

#label DE genes
et$DE = ""
et$DE[et$FDR < 0.05 & et$logFC > 0 ] = "Up"
et$DE[et$FDR < 0.05 & et$logFC < 0 ] = "Down"
et$DE[ et$DE == "" ] = "No"

table(et$DE)

write.table(et, "output/DiffExp_3xTgAD_BL6_hipp_F_18mon.txt", sep="\t", row.names = F, quote = F)

et$label = et$gene_name
for (i in c(1:dim(et)[1])) {
  if (et$DE[i] == "No") {
    et$label[i] = ""
  }
}


pdf("volcano plot/volcano_3xTgAD_BL6_hipp_F_18mon.pdf", width = 6, height = 6)
ggplot(data=et, aes(x=logFC, y=-log10(FDR), col=DE, label=label)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0, 0), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
dev.off()
