library("DOSE")
library(enrichplot)
library(org.Hs.eg.db)
library(org.Mmu.eg.db)
library(clusterProfiler)
#library(ReactomePA)
library(tidyverse)

rm(list =ls())

#import DEGs dataset
hspc_markers_240223 <- read_csv("./hsc_agetissue_markers_allgenes_240223.csv")
data <- hspc_markers_240223

#data cleaning
data<- data[data$logfoldchanges>0.5,]
data<- data[data$pvals_adj<0.05,]
data$group <- factor(data$group, levels = c('Early2nd_FL', 'Late2nd_FL', 'Early3rd_FL', 'Late2nd_BM', 'Early3rd_BM', 'Adult_BM'))
data <- data[order(data$group),]


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
write.csv(godf,file='240302_comparecluster_HSC_agetissut.csv')

#plot dotplot
fig <- dotplot(GO_ALL)  +
  theme(axis.text.y = element_text(size = 8, lineheight = 0.6),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 12))
fig

fig <- clusterProfiler::dotplot(GO_ALL,showCategory = 5,font.size  =10, label_format=100) +
  ggtitle("GO enrichment-ALL") +
  theme(axis.text.y = element_text(size = 10, lineheight = 0.6), axis.text.x = element_text(angle = 30, hjust = 1, size = 15)) +
  theme(legend.title = element_text(size=10), legend.text = element_text(size = 10), plot.title = element_text(size = 10) )
fig
ggsave(filename = "240308_compareCluster_GO.pdf",width = 8, height = 5,plot = fig)
