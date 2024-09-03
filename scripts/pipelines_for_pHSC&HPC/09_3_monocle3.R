# Trajectory analysis
## monocle3

library(monocle3)
library(ggplot2)
library(dplyr)

rm(list =ls())

# This is for MEMP pathway analysis
# Do the same thing for the other pathways

# import data
cds <- readRDS("./RDS/Pseudotime_fdg_subset_240201.RDS")
fdg_coordinates <- read.csv("./RDS/Pseudotime_fdg_subset_240201.csv", row.names = 1)
cds@int_colData@listData[["reducedDims"]][["UMAP"]] <- as.matrix(fdg_coordinates)

# PCA
set.seed(26)
cds <- preprocess_cds(cds, num_dim = 100, method="PCA", verbose=TRUE)

# Clustering
cds <- cluster_cells(cds, reduction_method = "UMAP", verbose=TRUE, resolution = 1e-07)
cds <- learn_graph(cds, learn_graph_control = list(ncenter=100))

# order cells
get_earliest_principal_node <- function(cds, cell_type="HSC/MPP"){
  cell_ids <- which(colData(cds)[, "cell_type"] == cell_type)
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

#plot
p1 <- plot_cells(cds, color_cells_by = "partition")
p2 <- plot_cells(cds, color_cells_by = "cell_type")
p3 = plot_cells(cds, color_cells_by = "cell_type", cell_size = 2, trajectory_graph_color="black",  label_cell_groups=FALSE, label_leaves=TRUE, label_branch_points=TRUE)
p4 = plot_cells(cds, color_cells_by = "pseudotime", cell_size = 1, trajectory_graph_color="black",
                label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE)#, graph_label_size=1.5)

ggsave(plot=p1,"./figures_R/240201_partition.pdf", device = "pdf")
ggsave(plot=p2,"./figures_R/240201_celltype.pdf", device = "pdf")
ggsave(plot=p3, "./figures_R/240201_trajectory.pdf", width=12, height=10)
ggsave(plot=p4, "./figures_R/240202_pseudotime.pdf", width=12, height=10)

#Finding genes that change as a function of pseudotime
cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)

#collect the trajectory-variable genes into modules
pr_deg_ids <- row.names(subset(cds_pr_test_res, q_value < 0.05))
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(10^seq(-7,-1)))
write.csv(gene_module_df, "./results/240202_test_MEMP_overall_degs_gene_modules.csv")

# combine the DEG info df (pr_deg_ids) with modules (gene_module_df)
pr_deg_ids <- subset(cds_pr_test_res, q_value < 0.05)
pr_deg_ids$id = rownames(pr_deg_ids)
gene_module_df_final <- as.data.frame(gene_module_df)
merged <- merge(pr_deg_ids, gene_module_df_final, by="id")
head(merged)
write.csv(merged, "./results/240202_MEMP_monocle3_overall_merged_deg_info.csv")

# plot the aggregate module scores
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)),cell_group=colData(cds)$cell_type)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
write.csv(as.matrix(agg_mat), "./results/240202_MEMP_monocle3_overall_agg_mat.csv")
p5 <- pheatmap::pheatmap(agg_mat,scale="row", clustering_method="ward.D2")
ggsave(plot =p5, "./figures_R/240202_MEMP_monocle3_gene_modules_heatmap.pdf", width=12, height=10)

## (option) change the order of modules
final_order = c("Module 7", "Module 10", "Module 3", "Module 4", "Module 9", "Module 13","Module 1","Module 6",
                "Module 11", "Module 12", "Module 2", "Module 5","Module 8")
celltype_order=c('HSC/MPP', 'HSC/CLP','Prolif-MPP1','Prolif-MPP2','Pre Pro-B','Pro B','pDC')
final_agg_mat <- agg_mat[final_order,]
final_agg_mat <- final_agg_mat[,celltype_order]
p6 <- pheatmap::pheatmap(final_agg_mat, scale="row", clustering_method="ward.D2", cluster_rows = FALSE, cluster_cols = FALSE)
ggsave(plot = p6, "./figures_R/240202_MEMP_monocle3_gene_modules_heatmap_order.pdf", width=10, height=10)

write.csv(as.matrix(final_agg_mat), "./results/240202_MEMP_monocle3_overall_final_agg_mat.csv")

# pass gene_module_df to plot_cells()
p7 <- plot_cells(cds,
           genes=gene_module_df,# %>% filter(module %in% c(3,6,8)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
ggsave(plot =p7, "./figures_R/240202_MEMP_module_celltype.pdf", width=10, height=10)

# save metadata
p = plot_cells(cds, color_cells_by = "pseudotime", cell_size = 1.5, trajectory_graph_color="black", label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=1.5)
point.data <- ggplot_build(p)[["plot"]][["data"]]
colnames(point.data)
point.data["cell_color"]
write.csv(point.data["cell_color"], "./results/240202_MEMP_monocle3_pst_metadata.csv")

# save cds
save_monocle_objects(cds=cds, directory_path='./RDS/240202_MEMP_objects', comment='Stored 2024-02-02.')
