# set directory
#############################################################################################
data_dir = "./data/human_pancreas"
out_dir = "./result/human_pancreas"

reactome_dir = "./data/gsea_cut_off/human_pancreas/human_pancreas_c2"
reactome_out_dir = "./result/gsea_cut_off/human_pancreas/human_pancreas_c2"

go_dir = "./data/gsea_cut_off/human_pancreas/human_pancreas_c5"
go_out_dir = "./result/gsea_cut_off/human_pancreas/human_pancreas_c5"

validation_dir = "./data/gsea_cut_off/human_pancreas/human_pancreas_c7"
validation_out_dir = "./result/gsea_cut_off/human_pancreas/human_pancreas_c7"
#############################################################################################



#prepare environment
#############################################################################################
if (!require(Seurat)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("Seurat")
}

setwd("E:/single_cell_experence/GSEA/gseapy_sourcecode/R_script")

library(Seurat)
library(ggplot2)
library(dplyr)
#############################################################################################




# origin run
#############################################################################################
pancreas.data <- Read10X(data.dir = data_dir)
print(data_dir)
n_feature <- 1
prject_name = "origin"

pancreas <- CreateSeuratObject(counts = pancreas.data, project = prject_name, min.cells = 3, min.features = n_feature)
pancreas


VlnPlot(pancreas, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)


nFeature_RNA_max <- 12000
nFeature_RNA_max <- as.numeric(nFeature_RNA_max)
scale_factor <- 30000
scale_factor <- as.numeric(scale_factor)



#normalizing the data

pancreas <- NormalizeData(pancreas, normalization.method = "LogNormalize", scale.factor = scale_factor)  

pancreas <- NormalizeData(pancreas)


pancreas <- FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = nFeature_RNA_max)    

t_n <- 5
t_n <- as.numeric(t_n)
top_n <- head(VariableFeatures(pancreas), t_n)


plot2_1 <- VariableFeaturePlot(pancreas)
plot2_2 <- LabelPoints(plot = plot2_1, points = top_n, repel = TRUE)
CombinePlots(plots = list(plot2_1, plot2_2))  



#scaling the data
all.genes <- rownames(pancreas)
pancreas <- ScaleData(pancreas, features = all.genes)

#run pca
pancreas <- RunPCA(pancreas, features = VariableFeatures(object = pancreas))    
print(pancreas[["pca"]], dims = 1:5, nfeatures = 5) 


VizDimLoadings(pancreas, dims = 1:2, reduction = "pca")




DimPlot(pancreas, reduction = "pca")




# cells_counts <- readline("please input the number of cells to plot heatmap!(tips:500/2700) :")
# cells_counts <- as.numeric(cells_counts)
# 
# dim_counts <- readline("please input the first n dim to plot heatmap :")
# dim_counts <- as.numeric(dim_counts)
# 
# DimHeatmap(pancreas, dims = 1:dim_counts, cells = cells_counts, balanced = TRUE)    

ElbowPlot(pancreas)

umap_pca <- 20
umap_pca <- as.numeric(umap_pca)


resolution_num = 0.8
resolution_num <- as.numeric(resolution_num)

pancreas <- FindNeighbors(pancreas, dims = 1:umap_pca)
pancreas <- FindClusters(pancreas, resolution = resolution_num)    

pancreas <- RunUMAP(pancreas, dims = 1:umap_pca)    


DimPlot(pancreas, reduction = "umap")    


pancreas.markers <- FindAllMarkers(pancreas, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)    
pancreas.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


for (ss in top_n){
  print(VlnPlot(pancreas, features = ss))
  tem <- readline("press any key to continue!")
}    
for (ss in top_n){
  print(FeaturePlot(pancreas, features = ss))
  tem <- readline("press any key to continue!")
  
}

top_n_show <- pancreas.markers %>% group_by(cluster) %>% top_n(n = t_n, wt = avg_log2FC)
DoHeatmap(pancreas, features = top_n_show$gene) + NoLegend()

top_n_each_cluster_show <- unique(top_n_show$gene)
for (ss in top_n_each_cluster_show){
  print(VlnPlot(pancreas, features = ss))
  tem <- readline("press any key to continue!")
}  

for (ss in top_n_each_cluster_show){
  print(FeaturePlot(pancreas, features = ss))
  tem <- readline("press any key to continue!")
}



# write top n variable for each cluster
top_n_each_cluster <- paste(out_dir,"top_n_variable_for_each_cluster.csv",sep = "/")
write.csv(top_n_show, top_n_each_cluster, row.names = TRUE)



# write cluster information!
cluster_all <- paste(out_dir,"cluster_all.csv",sep = "/")
write.csv(pancreas@meta.data$seurat_clusters, cluster_all, row.names = FALSE)

# write cluster marker information
cluster_marker_all <- paste(out_dir,"cluster_marker.csv", sep = "/")
write.csv(pancreas.markers,cluster_marker_all)



#replot
cell_types_fn <- paste(data_dir,"cell_type.tsv",sep = "/")
cell_types <- read.table(cell_types_fn, stringsAsFactors = FALSE)
cell_types <- cell_types[,1]
cell_types_l <-  as.factor(cell_types)
dim(cell_types_l)
cell_types_l


pancreas@meta.data$cell_types <- cell_types_l
DimPlot(pancreas, reduction = "umap", group.by = "cell_types", label = TRUE) 

DimPlot(pancreas, reduction = "umap", group.by = "seurat_clusters", label = TRUE) 

top_n_show_ct <- pancreas.markers %>% group_by("cell_types") %>% top_n(n = t_n, wt = avg_log2FC)
DoHeatmap(pancreas, features = top_n_show_ct$gene) + NoLegend()


#############################################################################################





rm(list=ls())


# function run reactome
#############################################################################################

# set directory
#############################################################################################
data_dir = "./data/human_pancreas"
out_dir = "./result/human_pancreas"

reactome_dir = "./data/gsea_cut_off/human_pancreas/human_pancreas_c2"
reactome_out_dir = "./result/gsea_cut_off/human_pancreas/human_pancreas_c2"

go_dir = "./data/gsea_cut_off/human_pancreas/human_pancreas_c5"
go_out_dir = "./result/gsea_cut_off/human_pancreas/human_pancreas_c5"

validation_dir = "./data/gsea_cut_off/human_pancreas/human_pancreas_c7"
validation_out_dir = "./result/gsea_cut_off/human_pancreas/human_pancreas_c7"
#############################################################################################



pancreas.data <- Read10X(data.dir = reactome_dir)
print(reactome_dir)
n_feature = 1
n_feature <- as.numeric(n_feature)
prject_name = "pancreas_data"

pancreas <- CreateSeuratObject(counts = pancreas.data, project = prject_name, min.cells = 3, min.features = n_feature)
pancreas


VlnPlot(pancreas, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

nFeature_RNA_max <- 1400

nFeature_RNA_max <- as.numeric(nFeature_RNA_max)
scale_factor <- 4000
scale_factor <- as.numeric(scale_factor)


#normalizing the data

pancreas <- NormalizeData(pancreas, normalization.method = "LogNormalize", scale.factor = scale_factor)    

pancreas <- NormalizeData(pancreas)


pancreas <- FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = nFeature_RNA_max)    

t_n <- 5
t_n <- as.numeric(t_n)
top_n <- head(VariableFeatures(pancreas), t_n)


plot2_1 <- VariableFeaturePlot(pancreas)
plot2_2 <- LabelPoints(plot = plot2_1, points = top_n, repel = TRUE)
CombinePlots(plots = list(plot2_1, plot2_2))  



#scaling the data
all.genes <- rownames(pancreas)
pancreas <- ScaleData(pancreas, features = all.genes)

#run pca
pancreas <- RunPCA(pancreas, features = VariableFeatures(object = pancreas))    
print(pancreas[["pca"]], dims = 1:5, nfeatures = 5) 

VizDimLoadings(pancreas, dims = 1:2, reduction = "pca")


DimPlot(pancreas, reduction = "pca")



# cells_counts <- readline("please input the number of cells to plot heatmap!(tips:500/2700) :")
# cells_counts <- as.numeric(cells_counts)
# 
# dim_counts <- readline("please input the first n dim to plot heatmap :")
# dim_counts <- as.numeric(dim_counts)
# 
# DimHeatmap(pancreas, dims = 1:dim_counts, cells = cells_counts, balanced = TRUE)    


ElbowPlot(pancreas)
umap_pca <- 20
umap_pca <- as.numeric(umap_pca)


resolution_num = 1.2
resolution_num <- as.numeric(resolution_num)

pancreas <- FindNeighbors(pancreas, dims = 1:umap_pca)
pancreas <- FindClusters(pancreas, resolution = resolution_num)    

pancreas <- RunUMAP(pancreas, dims = 1:umap_pca)    


DimPlot(pancreas, reduction = "umap")    


pancreas.markers <- FindAllMarkers(pancreas, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)    

pancreas.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)



for (ss in top_n){
  print(VlnPlot(pancreas, features = ss))
  tem <- readline("press any key to continue!")
}    
for (ss in top_n){
  print(FeaturePlot(pancreas, features = ss))
  tem <- readline("press any key to continue!")
  
}

top_n_show <- pancreas.markers %>% group_by(cluster) %>% top_n(n = t_n, wt = avg_log2FC)
DoHeatmap(pancreas, features = top_n_show$gene) + NoLegend()

top_n_each_cluster_show <- unique(top_n_show$gene)
for (ss in top_n_each_cluster_show){
  print(VlnPlot(pancreas, features = ss))
  tem <- readline("press any key to continue!")
}  

for (ss in top_n_each_cluster_show){
  print(FeaturePlot(pancreas, features = ss))
  tem <- readline("press any key to continue!")
}



# write top n variable for each cluster
top_n_each_cluster <- paste(reactome_out_dir,"reactome_top_n_variable_for_each_cluster.csv",sep = "/")
write.csv(top_n_show, top_n_each_cluster, row.names = TRUE)



# write cluster information
cluster_all <- paste(reactome_out_dir,"reactome_cluster_all.csv",sep = "/")
write.csv(pancreas@meta.data$seurat_clusters, cluster_all, row.names = FALSE)

# write cluster marker
cluster_marker_all <- paste(reactome_out_dir,"reactome_cluster_marker.csv", sep = "/")
write.csv(pancreas.markers,cluster_marker_all)



#replot
cell_types_fn <- paste(data_dir,"cell_type.tsv",sep = "/")
cell_types <- read.table(cell_types_fn, stringsAsFactors = FALSE)
cell_types <- cell_types[,1]
cell_types_l <-  as.factor(cell_types)
dim(cell_types_l)
cell_types_l


pancreas@meta.data$cell_types <- cell_types_l
DimPlot(pancreas, reduction = "umap", group.by = "cell_types", label = TRUE) 

DimPlot(pancreas, reduction = "umap", group.by = "seurat_clusters", label = TRUE) 

top_n_show_ct <- pancreas.markers %>% group_by('cell_types') %>% top_n(n = t_n, wt = avg_log2FC)




#heatmap
library(ggplot2)
data_heatmap <- pancreas@assays$RNA@scale.data
data_row_name <- pancreas@assays$RNA@data@Dimnames[[1]]
data_col_name <- pancreas@assays$RNA@data@Dimnames[[2]]
row.names(data_heatmap) <- data_row_name
colnames(data_heatmap) <- data_col_name

heatmap(data_heatmap,cexRow=0.5)


#############################################################################################
rm(list=ls())


# function run go
#############################################################################################

# set directory
#############################################################################################
data_dir = "./data/human_pancreas"
out_dir = "./result/human_pancreas"

reactome_dir = "./data/gsea_cut_off/human_pancreas/human_pancreas_c2"
reactome_out_dir = "./result/gsea_cut_off/human_pancreas/human_pancreas_c2"

go_dir = "./data/gsea_cut_off/human_pancreas/human_pancreas_c5"
go_out_dir = "./result/gsea_cut_off/human_pancreas/human_pancreas_c5"

validation_dir = "./data/gsea_cut_off/human_pancreas/human_pancreas_c7"
validation_out_dir = "./result/gsea_cut_off/human_pancreas/human_pancreas_c7"
#############################################################################################









pancreas.data <- Read10X(data.dir = go_dir)
print(go_dir)
n_feature = 1
n_feature <- as.numeric(n_feature)
prject_name = "GO"
pancreas <- CreateSeuratObject(counts = pancreas.data, project = prject_name, min.cells = 3, min.features = n_feature)
pancreas
VlnPlot(pancreas, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)



nFeature_RNA_max <- 7000
nFeature_RNA_max <- as.numeric(nFeature_RNA_max)
scale_factor <- 20000
scale_factor <- as.numeric(scale_factor)



#normalizing the data
pancreas <- NormalizeData(pancreas, normalization.method = "LogNormalize", scale.factor = scale_factor)   
pancreas <- NormalizeData(pancreas)


pancreas <- FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = nFeature_RNA_max)    

t_n <- 5
t_n <- as.numeric(t_n)
top_n <- head(VariableFeatures(pancreas), t_n)


#scaling the data
all.genes <- rownames(pancreas)
pancreas <- ScaleData(pancreas, features = all.genes)

#run pca
pancreas <- RunPCA(pancreas, features = VariableFeatures(object = pancreas))    
print(pancreas[["pca"]], dims = 1:5, nfeatures = 5) 


ElbowPlot(pancreas)
umap_pca <- 20
umap_pca <- as.numeric(umap_pca)


resolution_num = 0.8
resolution_num <- as.numeric(resolution_num)

pancreas <- FindNeighbors(pancreas, dims = 1:umap_pca)
pancreas <- FindClusters(pancreas, resolution = resolution_num)    

pancreas <- RunUMAP(pancreas, dims = 1:umap_pca)


DimPlot(pancreas, reduction = "umap")    


pancreas.markers <- FindAllMarkers(pancreas, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)    

pancreas.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


for (ss in top_n){
  print(VlnPlot(pancreas, features = ss))
  tem <- readline("press any key to continue!")
}    
for (ss in top_n){
  print(FeaturePlot(pancreas, features = ss))
  tem <- readline("press any key to continue!")
  
}

top_n_show <- pancreas.markers %>% group_by(cluster) %>% top_n(n = t_n, wt = avg_log2FC)
# DoHeatmap(pancreas, features = top_n_show$gene) + NoLegend()


top_n_each_cluster_show <- unique(top_n_show$gene)
for (ss in top_n_each_cluster_show){
  print(VlnPlot(pancreas, features = ss))
  tem <- readline("press any key to continue!")
}  

for (ss in top_n_each_cluster_show){
  print(FeaturePlot(pancreas, features = ss))
  tem <- readline("press any key to continue!")
  
}



#write top n variable for each cluster

top_n_each_cluster <- paste(go_out_dir,"go_top_n_variable_for_each_cluster.csv",sep = "/")
write.csv(top_n_show, top_n_each_cluster, row.names = TRUE)



#write cluster information

cluster_all <- paste(go_out_dir,"go_cluster_all.csv",sep = "/")

write.csv(pancreas@meta.data$seurat_clusters, cluster_all, row.names = FALSE)

# write marker information
cluster_marker_all <- paste(go_out_dir,"go_cluster_marker.csv", sep = "/")

write.csv(pancreas.markers,cluster_marker_all)


#replot
cell_types_fn <- paste(data_dir,"cell_type.tsv",sep = "/")
cell_types <- read.table(cell_types_fn, stringsAsFactors = FALSE)
cell_types <- cell_types[,1]
cell_types_l <-  as.factor(cell_types)
dim(cell_types_l)
cell_types_l


pancreas@meta.data$cell_types <- cell_types_l
DimPlot(pancreas, reduction = "umap", group.by = "cell_types", label = TRUE) 
DimPlot(pancreas, reduction = "umap", group.by = "cell_types", label = FALSE) 

DimPlot(pancreas, reduction = "umap", group.by = "seurat_clusters", label = TRUE) 
DimPlot(pancreas, reduction = "umap", group.by = "seurat_clusters", label = FALSE)

top_n_show_ct <- pancreas.markers %>% group_by("cell_types") %>% top_n(n = t_n, wt = avg_log2FC)



#heatmap
library(ggplot2)
data_heatmap <- pancreas@assays$RNA@scale.data
data_row_name <- pancreas@assays$RNA@data@Dimnames[[1]]
data_col_name <- pancreas@assays$RNA@data@Dimnames[[2]]
row.names(data_heatmap) <- data_row_name
colnames(data_heatmap) <- data_col_name

heatmap(data_heatmap,cexRow=0.5)

#############################################################################################





rm(list=ls())
# function on origin Reactome mm!
#############################################################################################

# set directory
#############################################################################################
data_dir = "./data/human_pancreas"
out_dir = "./result/human_pancreas"

reactome_dir = "./data/gsea_cut_off/human_pancreas/human_pancreas_c2"
reactome_out_dir = "./result/gsea_cut_off/human_pancreas/human_pancreas_c2"

go_dir = "./data/gsea_cut_off/human_pancreas/human_pancreas_c5"
go_out_dir = "./result/gsea_cut_off/human_pancreas/human_pancreas_c5"

validation_dir = "./data/gsea_cut_off/human_pancreas/human_pancreas_c7"
validation_out_dir = "./result/gsea_cut_off/human_pancreas/human_pancreas_c7"
#############################################################################################
s_data_dir = data_dir


pancreas.data <- Read10X(data.dir = s_data_dir)

n_feature = 1
n_feature <- as.numeric(n_feature)

prject_name = "REACTOME"
pancreas <- CreateSeuratObject(counts = pancreas.data, project = "ORINGIN", min.cells = 3, min.features = n_feature)

VlnPlot(pancreas, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

nFeature_RNA_max <- 12000

nFeature_RNA_max <- as.numeric(nFeature_RNA_max)

scale_factor <- 30000

scale_factor <- as.numeric(scale_factor)


#normalizing the data
pancreas <- NormalizeData(pancreas, normalization.method = "LogNormalize", scale.factor = scale_factor)   
pancreas <- NormalizeData(pancreas)  

pancreas <- FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = nFeature_RNA_max)    

t_n <- 5
t_n <- as.numeric(t_n)
top_n <- head(VariableFeatures(pancreas), t_n)

# plot2_1 <- VariableFeaturePlot(pancreas)
# plot2_2 <- LabelPoints(plot = plot2_1, points = top_n, repel = TRUE)
# CombinePlots(plots = list(plot2_1, plot2_2))     



#scaling the data
all.genes <- rownames(pancreas)
pancreas <- ScaleData(pancreas, features = all.genes)

#run pca
pancreas <- RunPCA(pancreas, features = VariableFeatures(object = pancreas))    
print(pancreas[["pca"]], dims = 1:5, nfeatures = 5)


VizDimLoadings(pancreas, dims = 1:2, reduction = "pca")


DimPlot(pancreas, reduction = "pca")


ElbowPlot(pancreas)

umap_pca <- 20
umap_pca <- as.numeric(umap_pca)

resolution_num = 0.8
resolution_num <- as.numeric(resolution_num)

pancreas <- FindNeighbors(pancreas, dims = 1:umap_pca)
pancreas <- FindClusters(pancreas, resolution = resolution_num)   

pancreas <- RunUMAP(pancreas, dims = 1:umap_pca)   

DimPlot(pancreas, reduction = "umap")    


pancreas.markers <- FindAllMarkers(pancreas, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pancreas.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

DimPlot(pancreas, label = TRUE) + NoLegend()



#Add a function Seurat Object to plot!
pancreas_func.data <- Read10X(data.dir = reactome_dir)

n_feature_func = 2
n_feature_func <- as.numeric(n_feature_func)
func_obj_name = "function name!"
func_obj <- CreateSeuratObject(counts = pancreas_func.data, project = func_obj_name, min.cells = 3, min.features = n_feature_func)


VlnPlot(func_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
pancreas[["function_n"]] = CreateAssayObject(pancreas_func.data)


func_feature_max <- 1500

func_feature_max <- as.numeric(func_feature_max)

func_scale_factor <- 4000
func_scale_factor <- as.numeric(func_scale_factor)

pancreas <- NormalizeData(pancreas, assay = "function_n", normalization.method = "LogNormalize",  scale.factor = scale_factor)

pancreas <- ScaleData(pancreas, assay = "function_n")

pancreas <- FindVariableFeatures(pancreas, assay = "function_n", selection.method = "vst", nfeatures = func_feature_max)

t_n_f <- 20
t_n_f <- as.numeric(t_n_f)
top_n_f <- head(VariableFeatures(pancreas, assay = "function_n"), t_n_f)

special_function <- "reactome-response-to-elevated-platelet-cytosolic-ca2plus"
print(FeaturePlot(pancreas, features = special_function))
print(RidgePlot(pancreas, features = special_function))

for (ss in top_n_f){
  print(FeaturePlot(pancreas, features = ss))
  tem <- readline("press any key to continue!")
}

for (ss in top_n_f){
  print(RidgePlot(pancreas, features = ss))
  tem <- readline("press any key to continue!")
}


# #Identify differentially expressed functions between clusters
# 
# max_cells = 777
# 
# pancreas.small <- subset(pancreas, downsample = max_cells)
# 
# functions.markers <- FindAllMarkers(pancreas , assay = "function_n", only_pos = TRUE)
# 
# DoHeatmap(pancreas.small, features = unique(functions.markers$gene), assay = "function_n", angle = 90) + NoLegend()


pancreas_func.markers <- FindAllMarkers(pancreas, assay = "function_n", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pancreas_func.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> res_cluster_top_gene



top_n_function_show <- pancreas_func.markers %>% group_by(cluster) %>% top_n(n = t_n, wt = avg_log2FC)
#write top n variable for each cluster

top_n_each_cluster <- paste(reactome_out_dir,"reactome_mm_top_n_variable_funciotn_for_each_cluster.csv",sep = "/")
write.csv(top_n_function_show, top_n_each_cluster, row.names = TRUE)


cluster_marker_gene <- unique(res_cluster_top_gene$gene)

for (ss in cluster_marker_gene){
  print(FeaturePlot(pancreas, features = ss, label = TRUE))
  tem <- readline("press any key to continue!")
}

for (ss in cluster_marker_gene){
  print(RidgePlot(pancreas, features = ss))
  tem <- readline("press any key to continue!")
}


# write cluster information!
cluster_all <- paste(reactome_out_dir, "reactome_mm_cluster_all.csv", sep = "/")
write.csv(pancreas@meta.data$seurat_clusters, cluster_all, row.names = FALSE)

#write cluster marker (gene)
cluster_marker_source <- paste(reactome_out_dir,"reactome_mm_cluster_marker_source.csv", sep = "/")
write.csv(pancreas.markers,cluster_marker_source)

#write cluster marker (function)
cluster_marker_function <- paste(reactome_out_dir, "reactome_mm_cluster_marker_function.csv", sep = "/")
write.csv(pancreas_func.markers,cluster_marker_function)


#replot
cell_types_fn <- paste(data_dir,"cell_type.tsv",sep = "/")
cell_types <- read.table(cell_types_fn, stringsAsFactors = FALSE)
cell_types <- cell_types[,1]
cell_types_l <-  as.factor(cell_types)
dim(cell_types_l)
cell_types_l


pancreas@meta.data$cell_types <- cell_types_l
DimPlot(pancreas, reduction = "umap", group.by = "cell_types", label = TRUE) 

DimPlot(pancreas, reduction = "umap", group.by = "seurat_clusters", label = TRUE) 

top_n_show_ct <- pancreas.markers %>% group_by("cell_types") %>% top_n(n = t_n, wt = avg_log2FC)
DoHeatmap(pancreas, features = top_n_show_ct$gene) + NoLegend()



#############################################################################################
rm(list=ls())



# function on origin go mm!
#############################################################################################
# set directory
#############################################################################################
data_dir = "./data/human_pancreas"
out_dir = "./result/human_pancreas"

reactome_dir = "./data/gsea_cut_off/human_pancreas/human_pancreas_c2"
reactome_out_dir = "./result/gsea_cut_off/human_pancreas/human_pancreas_c2"

go_dir = "./data/gsea_cut_off/human_pancreas/human_pancreas_c5"
go_out_dir = "./result/gsea_cut_off/human_pancreas/human_pancreas_c5"

validation_dir = "./data/gsea_cut_off/human_pancreas/human_pancreas_c7"
validation_out_dir = "./result/gsea_cut_off/human_pancreas/human_pancreas_c7"
#############################################################################################




s_data_dir <- data_dir


pancreas.data <- Read10X(data.dir = s_data_dir)

n_feature = 1
n_feature <- as.numeric(n_feature)

prject_name = "go"
pancreas <- CreateSeuratObject(counts = pancreas.data, project = "ORINGIN", min.cells = 3, min.features = n_feature)
pancreas
VlnPlot(pancreas, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)


nFeature_RNA_max <- 12000
nFeature_RNA_max <- as.numeric(nFeature_RNA_max)
scale_factor <- 30000
scale_factor <- as.numeric(scale_factor)


#normalizing the data
pancreas <- NormalizeData(pancreas, normalization.method = "LogNormalize", scale.factor = scale_factor)   
pancreas <- NormalizeData(pancreas)  

pancreas <- FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = nFeature_RNA_max)    

t_n <- 5
t_n <- as.numeric(t_n)
top_n <- head(VariableFeatures(pancreas), t_n)

# plot2_1 <- VariableFeaturePlot(pancreas)
# plot2_2 <- LabelPoints(plot = plot2_1, points = top_n, repel = TRUE)
# CombinePlots(plots = list(plot2_1, plot2_2))     


#scaling the data
all.genes <- rownames(pancreas)
pancreas <- ScaleData(pancreas, features = all.genes)

#run pca
pancreas <- RunPCA(pancreas, features = VariableFeatures(object = pancreas))    
print(pancreas[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pancreas, dims = 1:2, reduction = "pca")


DimPlot(pancreas, reduction = "pca")


ElbowPlot(pancreas)
umap_pca <- 20
umap_pca <- as.numeric(umap_pca)

resolution_num = 0.8
resolution_num <- as.numeric(resolution_num)

pancreas <- FindNeighbors(pancreas, dims = 1:umap_pca)
pancreas <- FindClusters(pancreas, resolution = resolution_num)   

pancreas <- RunUMAP(pancreas, dims = 1:umap_pca)   

DimPlot(pancreas, reduction = "umap")


pancreas.markers <- FindAllMarkers(pancreas, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pancreas.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

DimPlot(pancreas, label = TRUE) + NoLegend()



#Add a function Seurat Object to plot!
pancreas_func.data <- Read10X(data.dir = go_dir)

n_feature_func = 2
n_feature_func <- as.numeric(n_feature_func)
func_obj_name = "function name!"
func_obj <- CreateSeuratObject(counts = pancreas_func.data, project = func_obj_name, min.cells = 3, min.features = n_feature_func)



VlnPlot(func_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)


pancreas[["function_n"]] = CreateAssayObject(pancreas_func.data)


func_feature_max <- 7000

func_feature_max <- as.numeric(func_feature_max)

func_scale_factor <- 20000
func_scale_factor <- as.numeric(func_scale_factor)

pancreas <- NormalizeData(pancreas, assay = "function_n", normalization.method = "LogNormalize",  scale.factor = scale_factor)

pancreas <- ScaleData(pancreas, assay = "function_n")

pancreas <- FindVariableFeatures(pancreas, assay = "function_n", selection.method = "vst", nfeatures = func_feature_max)

t_n_f <- 20
t_n_f <- as.numeric(t_n_f)
top_n_f <- head(VariableFeatures(pancreas, assay = "function_n"), t_n_f)


for (ss in top_n_f){
  print(FeaturePlot(pancreas, features = ss))
  tem <- readline("press any key to continue!")
}

for (ss in top_n_f){
  print(RidgePlot(pancreas, features = ss))
  tem <- readline("press any key to continue!")
}


#Identify differentially expressed functions between clusters

# max_cells = 2126
# 
# pancreas.small <- subset(pancreas, downsample = max_cells)
# 
# functions.markers <- FindAllMarkers(pancreas , assay = "function_n", only_pos = TRUE)
# 
# DoHeatmap(pancreas.small, features = unique(functions.markers$gene), assay = "function_n", angle = 90) + NoLegend()
# 
# tem <- readline("press any key to continue!") 

pancreas_func.markers <- FindAllMarkers(pancreas, assay = "function_n", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pancreas_func.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> res_cluster_top_gene

DoHeatmap(pancreas.small, features = unique(res_cluster_top_gene$gene), assay = "function_n", angle = 90) + NoLegend()


top_n_function_show <- pancreas_func.markers %>% group_by(cluster) %>% top_n(n = t_n, wt = avg_log2FC)
#write top n variable for each cluster

top_n_each_cluster <- paste(go_out_dir,"go_mm_top_n_variable_funciotn_for_each_cluster.csv",sep = "/")
write.csv(top_n_function_show, top_n_each_cluster, row.names = TRUE)







cluster_marker_gene <- unique(res_cluster_top_gene$gene)

for (ss in cluster_marker_gene){
  print(FeaturePlot(pancreas, features = ss, label = TRUE))
  dev.copy(which=dev_n)
  dev.set(which=dev_default)
  tem <- readline("press any key to continue!")
}

for (ss in cluster_marker_gene){
  print(RidgePlot(pancreas, features = ss))
  dev.copy(which=dev_n)
  dev.set(which=dev_default)
  tem <- readline("press any key to continue!")
}



#write some information to output file!
cluster_all <- paste(go_out_dir, "go_mm_cluster_all.csv", sep = "/")

write.csv(pancreas@meta.data$seurat_clusters, cluster_all, row.names = FALSE)


cluster_marker_source <- paste(go_out_dir,"go_mm_cluster_marker_source.csv", sep = "/")

write.csv(pancreas.markers,cluster_marker_source)

cluster_marker_function <- paste(go_out_dir, "go_mm_cluster_marker_function.csv", sep = "/")

write.csv(pancreas_func.markers,cluster_marker_function)



#replot
cell_types_fn <- paste(data_dir,"cell_type.tsv",sep = "/")
cell_types <- read.table(cell_types_fn, stringsAsFactors = FALSE)
cell_types <- cell_types[,1]
cell_types_l <-  as.factor(cell_types)
dim(cell_types_l)
cell_types_l


pancreas@meta.data$cell_types <- cell_types_l
DimPlot(pancreas, reduction = "umap", group.by = "cell_types", label = TRUE) 
DimPlot(pancreas, reduction = "umap", group.by = "cell_types", label = FALSE) 


DimPlot(pancreas, reduction = "umap", group.by = "seurat_clusters", label = TRUE) 
DimPlot(pancreas, reduction = "umap", group.by = "seurat_clusters", label = FALSE) 



top_n_show_ct <- pancreas_func.markers %>% group_by("cell_types") %>% top_n(n = t_n, wt = avg_log2FC)
DoHeatmap(pancreas, features = top_n_show_ct$gene) + NoLegend()
dev.copy(which=dev_n)
dev.set(which=dev_default)




rm(list=ls())
dev.off()
dev.off()

#############################################################################################