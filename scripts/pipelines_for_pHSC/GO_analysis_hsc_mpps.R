library("DOSE")
library(enrichplot)
library(org.Hs.eg.db)
library(org.Mmu.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(tidyverse)

rm(list =ls())

#import DEGs dataset
hspc_markers_240130 <- read_csv("~/Documents/Mf_10x/240109_ALL_HSC/results/hspc_markers_allgenes_240130.csv")
data <- hspc_markers_240130

#data cleaning
data<- data[data$logfoldchanges>0.5,]
data<- data[data$pvals_adj<0.05,]
data<- data[data$pct_nz_group>0.01,]


#Gene name conversion
gene_list <- split(data$names,data$group)

for(i in names(gene_list)){
  gene_list[[i]] <- bitr(gene_list[[i]],
                         fromType = "SYMBOL",
                         toType = "ENTREZID",
                         OrgDb = "org.Mmu.eg.db")
}

for(i in names(gene_list)){
  gene_list[[i]] <- gene_list[[i]]$ENTREZID
}

#GO analysis
GO_ALL <- compareCluster(geneClusters = gene_list,
                         fun = "enrichGO",
                         ont="ALL",
                         OrgDb = "org.Mmu.eg.db")
godf <- as.data.frame(GO_ALL@compareClusterResult)
write.csv(godf,file='240302_comparecluster_HSC_MPP.csv')

#plot dotplot
fig <- dotplot(GO_ALL)  +
  theme(axis.text.y = element_text(size = 8, lineheight = 0.6),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 12))
fig
ggsave(filename = "240130_compareCluster_GO_ALL_logFC0.5.pdf",width = 8, height = 6,plot = fig)

fig <- clusterProfiler::dotplot(GO_ALL,showCategory = 10,font.size  =14, label_format=40) + ggtitle("GO enrichment-ALL") +
  theme(legend.title = element_text(size=10), legend.text = element_text(size = 10), plot.title = element_text(size = 10) )
fig
