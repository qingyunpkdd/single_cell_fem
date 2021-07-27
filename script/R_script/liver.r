# set directory
#############################################################################################
data_dir = "./data/human_liver"
out_dir = "./result/human_liver"

reactome_dir = "./data/gsea_cut_off/human_liver/human_liver_c2"
reactome_out_dir = "./result/gsea_cut_off/human_liver/human_liver_c2"

go_dir = "./data/gsea_cut_off/human_liver/human_liver_c5"
go_out_dir = "./result/gsea_cut_off/human_liver/human_liver_c5"

validation_dir = "./data/gsea_cut_off/human_liver/human_liver_c7"
validation_out_dir = "./result/gsea_cut_off/human_liver/human_liver_c7"
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
liver.data <- Read10X(data.dir = data_dir)
print(data_dir)
n_feature <- 1
prject_name = "origin"

liver <- CreateSeuratObject(counts = liver.data, project = prject_name, min.cells = 3, min.features = n_feature)
liver

VlnPlot(liver, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
nFeature_RNA_max <- 8000
nFeature_RNA_max <- as.numeric(nFeature_RNA_max)
scale_factor <- 35000
scale_factor <- as.numeric(scale_factor)



#normalizing the data
liver <- NormalizeData(liver, normalization.method = "LogNormalize", scale.factor = scale_factor)  
liver <- NormalizeData(liver)


liver <- FindVariableFeatures(liver, selection.method = "vst", nfeatures = nFeature_RNA_max)    
t_n <- 10
t_n <- as.numeric(t_n)
top_n <- head(VariableFeatures(liver), t_n)

# tem <- readline("press any key to continue!")   

#scaling the data
all.genes <- rownames(liver)
liver <- ScaleData(liver, features = all.genes)

#run pca
liver <- RunPCA(liver, features = VariableFeatures(object = liver))    
print(liver[["pca"]], dims = 1:5, nfeatures = 5)

umap_pca <- 15
umap_pca <- as.numeric(umap_pca)


resolution_num = 0.5
resolution_num <- as.numeric(resolution_num)

liver <- FindNeighbors(liver, dims = 1:umap_pca)
liver <- FindClusters(liver, resolution = resolution_num)    

liver <- RunUMAP(liver, dims = 1:umap_pca)    
liver.markers <- FindAllMarkers(liver, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)    
liver.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


for (ss in top_n){
  print(VlnPlot(liver, features = ss))
  tem <- readline("press any key to continue!")
}    
for (ss in top_n){
  print(FeaturePlot(liver, features = ss))
  tem <- readline("press any key to continue!")
  
}

top_n_show <- liver.markers %>% group_by(cluster) %>% top_n(n = t_n, wt = avg_log2FC)
DoHeatmap(liver, features = top_n_show$gene) + NoLegend()
tem <- readline("press any key to continue!")

top_n_each_cluster_show <- unique(top_n_show$gene)
for (ss in top_n_each_cluster_show){
  print(VlnPlot(liver, features = ss))
  tem <- readline("press any key to continue!")
}  

for (ss in top_n_each_cluster_show){
  print(FeaturePlot(liver, features = ss))
  tem <- readline("press any key to continue!")
  
}



#write top n variable for each cluster
top_n_each_cluster <- paste(out_dir,"top_n_variable_for_each_cluster.csv",sep = "/")
write.csv(top_n_show, top_n_each_cluster, row.names = TRUE)


# write cluster information for each clusters
cluster_all <- paste(out_dir,"cluster_all.csv",sep = "/")
write.csv(liver@meta.data$seurat_clusters, cluster_all, row.names = FALSE)


# write cluster markers
cluster_marker_all <- paste(out_dir,"cluster_marker.csv", sep = "/")
write.csv(liver.markers,cluster_marker_all)



#replot
cell_types_fn <- paste(data_dir,"cell_type.tsv",sep = "/")
cell_types <- read.table(cell_types_fn, stringsAsFactors = FALSE)
cell_types <- cell_types[,1]
cell_types_l <-  as.factor(cell_types)
dim(cell_types_l)
cell_types_l


liver@meta.data$cell_types <- cell_types_l
DimPlot(liver, reduction = "umap", group.by = "cell_types", label = TRUE) 

DimPlot(liver, reduction = "umap", group.by = "seurat_clusters", label = TRUE) 

top_n_show_ct <- liver.markers %>% group_by("cell_types") %>% top_n(n = t_n, wt = avg_log2FC)
DoHeatmap(liver, features = top_n_show_ct$gene) + NoLegend()



#############################################################################################


rm(list=ls())




# function run reactome
#############################################################################################
# set directory
#############################################################################################
data_dir = "./data/human_liver"
out_dir = "./result/human_liver"

reactome_dir = "./data/gsea_cut_off/human_liver/human_liver_c2"
reactome_out_dir = "./result/gsea_cut_off/human_liver/human_liver_c2"

go_dir = "./data/gsea_cut_off/human_liver/human_liver_c5"
go_out_dir = "./result/gsea_cut_off/human_liver/human_liver_c5"

validation_dir = "./data/gsea_cut_off/human_liver/human_liver_c7"
validation_out_dir = "./result/gsea_cut_off/human_liver/human_liver_c7"
#############################################################################################

liver.data <- Read10X(data.dir = reactome_dir)
print(reactome_dir)
n_feature = 1
n_feature <- as.numeric(n_feature)
prject_name = "reactome"

liver <- CreateSeuratObject(counts = liver.data, project = prject_name, min.cells = 3, min.features = n_feature)
liver
VlnPlot(liver, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

nFeature_RNA_max <- 1400
nFeature_RNA_max <- as.numeric(nFeature_RNA_max)
scale_factor <- 4500
scale_factor <- as.numeric(scale_factor)


#normalizing the data
liver <- NormalizeData(liver, normalization.method = "LogNormalize", scale.factor = scale_factor)    
liver <- NormalizeData(liver)


liver <- FindVariableFeatures(liver, selection.method = "vst", nfeatures = nFeature_RNA_max)    
t_n <- readline("please input the top n variable feature:")
t_n <- as.numeric(t_n)
top_n <- head(VariableFeatures(liver), t_n)



#scaling the data
all.genes <- rownames(liver)
liver <- ScaleData(liver, features = all.genes)

#run pca
liver <- RunPCA(liver, features = VariableFeatures(object = liver))    
print(liver[["pca"]], dims = 1:5, nfeatures = 5) 


umap_pca <- 20 
umap_pca <- as.numeric(umap_pca)


resolution_num = 0.6
resolution_num <- as.numeric(resolution_num)

liver <- FindNeighbors(liver, dims = 1:umap_pca)
liver <- FindClusters(liver, resolution = resolution_num)    

liver <- RunUMAP(liver, dims = 1:umap_pca)    


DimPlot(liver, reduction = "umap")    

liver.markers <- FindAllMarkers(liver, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)    

liver.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)






for (ss in top_n){
  print(VlnPlot(liver, features = ss))
  tem <- readline("press any key to continue!")
}    
for (ss in top_n){
  print(FeaturePlot(liver, features = ss))
  tem <- readline("press any key to continue!")
  
}

top_n_show <- liver.markers %>% group_by(cluster) %>% top_n(n = t_n, wt = avg_log2FC)
DoHeatmap(liver, features = top_n_show$gene) + NoLegend()


top_n_each_cluster_show <- unique(top_n_show$gene)
for (ss in top_n_each_cluster_show){
  print(VlnPlot(liver, features = ss))
  tem <- readline("press any key to continue!")
}  

for (ss in top_n_each_cluster_show){
  print(FeaturePlot(liver, features = ss))
  tem <- readline("press any key to continue!")
}



#write top n variable for each cluster

top_n_each_cluster <- paste(reactome_out_dir,"reactome_top_n_variable_for_each_cluster.csv",sep = "/")
write.csv(top_n_show, top_n_each_cluster, row.names = TRUE)



# write cluster information for each clusters

cluster_all <- paste(reactome_out_dir,"reactome_cluster_all.csv",sep = "/")

write.csv(liver@meta.data$seurat_clusters, cluster_all, row.names = FALSE)

#write cluster markers
cluster_marker_all <- paste(reactome_out_dir,"reactome_cluster_marker.csv", sep = "/")

write.csv(liver.markers,cluster_marker_all)



#replot
cell_types_fn <- paste(data_dir,"cell_type.tsv",sep = "/")
cell_types <- read.table(cell_types_fn, stringsAsFactors = FALSE)
cell_types <- cell_types[,1]
cell_types_l <-  as.factor(cell_types)
dim(cell_types_l)
cell_types_l


liver@meta.data$cell_types <- cell_types_l
DimPlot(liver, reduction = "umap", group.by = "cell_types", label = TRUE) 

DimPlot(liver, reduction = "umap", group.by = "seurat_clusters", label = TRUE) 

top_n_show_ct <- liver.markers %>% group_by(cell_types) %>% top_n(n = t_n, wt = avg_log2FC)
DoHeatmap(pbmc, features = top_n_show_ct$gene) + NoLegend()




#heatmap
library(ggplot2)
data_heatmap <- liver@assays$RNA@scale.data
data_row_name <- liver@assays$RNA@data@Dimnames[[1]]
data_col_name <- liver@assays$RNA@data@Dimnames[[2]]
row.names(data_heatmap) <- data_row_name
colnames(data_heatmap) <- data_col_name

heatmap(data_heatmap,cexRow=0.5)


#############################################################################################

rm(list=ls())


# function run go
#############################################################################################
# set directory
#############################################################################################
data_dir = "./data/human_liver"
out_dir = "./result/human_liver"

reactome_dir = "./data/gsea_cut_off/human_liver/human_liver_c2"
reactome_out_dir = "./result/gsea_cut_off/human_liver/human_liver_c2"

go_dir = "./data/gsea_cut_off/human_liver/human_liver_c5"
go_out_dir = "./result/gsea_cut_off/human_liver/human_liver_c5"

validation_dir = "./data/gsea_cut_off/human_liver/human_liver_c7"
validation_out_dir = "./result/gsea_cut_off/human_liver/human_liver_c7"
#############################################################################################





liver.data <- Read10X(data.dir = go_dir)
print(go_dir)
n_feature = 1
n_feature <- as.numeric(n_feature)
prject_name = "GO"
liver <- CreateSeuratObject(counts = liver.data, project = prject_name, min.cells = 3, min.features = n_feature)
liver
VlnPlot(liver, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

nFeature_RNA_max <- 6000
nFeature_RNA_max <- as.numeric(nFeature_RNA_max)
scale_factor <- 20000
scale_factor <- as.numeric(scale_factor)



#normalizing the data
liver <- NormalizeData(liver, normalization.method = "LogNormalize", scale.factor = scale_factor)    
liver <- NormalizeData(liver)

liver <- FindVariableFeatures(liver, selection.method = "vst", nfeatures = nFeature_RNA_max)    
t_n <- 5
t_n <- as.numeric(t_n)
top_n <- head(VariableFeatures(liver), t_n)


#scaling the data
all.genes <- rownames(liver)
liver <- ScaleData(liver, features = all.genes)

#run pca
liver <- RunPCA(liver, features = VariableFeatures(object = liver))    
print(liver[["pca"]], dims = 1:5, nfeatures = 5) 


umap_pca <- 15
umap_pca <- as.numeric(umap_pca)

# ?
resolution_num = 0.8
resolution_num <- as.numeric(resolution_num)

liver <- FindNeighbors(liver, dims = 1:umap_pca)
liver <- FindClusters(liver, resolution = resolution_num)    
liver <- RunUMAP(liver, dims = 1:umap_pca)    


# DimPlot(liver, reduction = "umap")    
 
liver.markers <- FindAllMarkers(liver, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)    

liver.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


for (ss in top_n){
  print(VlnPlot(liver, features = ss))
  tem <- readline("press any key to continue!")
}    
for (ss in top_n){
  print(FeaturePlot(liver, features = ss))
  tem <- readline("press any key to continue!")
  
}

top_n_show <- liver.markers %>% group_by(cluster) %>% top_n(n = t_n, wt = avg_log2FC)
# DoHeatmap(liver, features = top_n_show$gene) + NoLegend()

top_n_each_cluster_show <- unique(top_n_show$gene)
for (ss in top_n_each_cluster_show){
  print(VlnPlot(liver, features = ss))
  tem <- readline("press any key to continue!")
}  

for (ss in top_n_each_cluster_show){
  print(FeaturePlot(liver, features = ss))
  tem <- readline("press any key to continue!")
  
}



#write top n variable for each cluster
top_n_each_cluster <- paste(go_out_dir,"go_top_n_variable_for_each_cluster.csv",sep = "/")
write.csv(top_n_show, top_n_each_cluster, row.names = TRUE)



# write cluster information for each clusters
cluster_all <- paste(go_out_dir,"go_cluster_all.csv",sep = "/")
write.csv(liver@meta.data$seurat_clusters, cluster_all, row.names = FALSE)


#write cluster markers
cluster_marker_all <- paste(go_out_dir,"go_cluster_marker.csv", sep = "/")
write.csv(liver.markers,cluster_marker_all)


#replot
cell_types_fn <- paste(data_dir,"cell_type.tsv",sep = "/")
cell_types <- read.table(cell_types_fn, stringsAsFactors = FALSE)
cell_types <- cell_types[,1]
cell_types_l <-  as.factor(cell_types)
dim(cell_types_l)
cell_types_l


liver@meta.data$cell_types <- cell_types_l
DimPlot(liver, reduction = "umap", group.by = "cell_types", label = TRUE) 
DimPlot(liver, reduction = "umap", group.by = "cell_types", label = FALSE) 

DimPlot(liver, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(liver, reduction = "umap", group.by = "seurat_clusters", label = FALSE) 

top_n_show_ct <- liver.markers %>% group_by("cell_types") %>% top_n(n = t_n, wt = avg_log2FC)
DoHeatmap(liver, features = top_n_show_ct$gene) + NoLegend()



#heatmap
library(ggplot2)
data_heatmap <- liver@assays$RNA@scale.data
data_row_name <- liver@assays$RNA@data@Dimnames[[1]]
data_col_name <- liver@assays$RNA@data@Dimnames[[2]]
row.names(data_heatmap) <- data_row_name
colnames(data_heatmap) <- data_col_name

heatmap(data_heatmap,cexRow=0.5)

#############################################################################################


rm(list=ls())


# function on origin Reactome mm !
#############################################################################################
# set directory
#############################################################################################
data_dir = "./data/human_liver"
out_dir = "./result/human_liver"

reactome_dir = "./data/gsea_cut_off/human_liver/human_liver_c2"
reactome_out_dir = "./result/gsea_cut_off/human_liver/human_liver_c2"

go_dir = "./data/gsea_cut_off/human_liver/human_liver_c5"
go_out_dir = "./result/gsea_cut_off/human_liver/human_liver_c5"

validation_dir = "./data/gsea_cut_off/human_liver/human_liver_c7"
validation_out_dir = "./result/gsea_cut_off/human_liver/human_liver_c7"
#############################################################################################

s_data_dir <- data_dir


liver.data <- Read10X(data.dir = s_data_dir)

n_feature = 1
n_feature <- as.numeric(n_feature)

prject_name = "REACTOME"
liver <- CreateSeuratObject(counts = liver.data, project = "ORINGIN", min.cells = 3, min.features = n_feature)

VlnPlot(liver, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)



nFeature_RNA_max <- 8000
nFeature_RNA_max <- as.numeric(nFeature_RNA_max)
scale_factor <- 35000
scale_factor <- as.numeric(scale_factor)


#normalizing the data
liver <- NormalizeData(liver, normalization.method = "LogNormalize", scale.factor = scale_factor)   
liver <- NormalizeData(liver)  

liver <- FindVariableFeatures(liver, selection.method = "vst", nfeatures = nFeature_RNA_max)    

t_n <- 5
t_n <- as.numeric(t_n)
top_n <- head(VariableFeatures(liver), t_n)

# plot2_1 <- VariableFeaturePlot(liver)
# plot2_2 <- LabelPoints(plot = plot2_1, points = top_n, repel = TRUE)
# CombinePlots(plots = list(plot2_1, plot2_2))     


#scaling the data
all.genes <- rownames(liver)
liver <- ScaleData(liver, features = all.genes)

#run pca
liver <- RunPCA(liver, features = VariableFeatures(object = liver))    
print(liver[["pca"]], dims = 1:5, nfeatures = 5)


VizDimLoadings(liver, dims = 1:2, reduction = "pca")

DimPlot(liver, reduction = "pca")


ElbowPlot(liver)
umap_pca <- 15
umap_pca <- as.numeric(umap_pca)

resolution_num =  0.5
resolution_num <- as.numeric(resolution_num)

liver <- FindNeighbors(liver, dims = 1:umap_pca)
liver <- FindClusters(liver, resolution = resolution_num)   

liver <- RunUMAP(liver, dims = 1:umap_pca)   

DimPlot(liver, reduction = "umap")    


liver.markers <- FindAllMarkers(liver, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
liver.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

DimPlot(liver, label = TRUE) + NoLegend()


#Add a function Seurat Object to plot!
liver_func.data <- Read10X(data.dir = reactome_dir)

n_feature_func = 2
n_feature_func <- as.numeric(n_feature_func)
func_obj_name = "function name!"
func_obj <- CreateSeuratObject(counts = liver_func.data, project = func_obj_name, min.cells = 3, min.features = n_feature_func)



VlnPlot(func_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)


liver[["function_n"]] = CreateAssayObject(liver_func.data)


func_feature_max <- 1400

func_feature_max <- as.numeric(func_feature_max)

func_scale_factor <- 4500
func_scale_factor <- as.numeric(func_scale_factor)

liver <- NormalizeData(liver, assay = "function_n", normalization.method = "LogNormalize",  scale.factor = scale_factor)

liver <- ScaleData(liver, assay = "function_n")

liver <- FindVariableFeatures(liver, assay = "function_n", selection.method = "vst", nfeatures = func_feature_max)

t_n_f <- 20
t_n_f <- as.numeric(t_n_f)
top_n_f <- head(VariableFeatures(liver, assay = "function_n"), t_n_f)

# special_function <- "reactome-response-to-elevated-platelet-cytosolic-ca2plus"
# print(FeaturePlot(liver, features = special_function))
# print(RidgePlot(liver, features = special_function))

for (ss in top_n_f){
  print(FeaturePlot(liver, features = ss))
  tem <- readline("press any key to continue!")
}

for (ss in top_n_f){
  print(RidgePlot(liver, features = ss))
  tem <- readline("press any key to continue!")
}


#Identify differentially expressed functions between clusters

# max_cells = 777
# 
# liver.small <- subset(liver, downsample = max_cells)
# 
functions.markers <- FindAllMarkers(liver , assay = "function_n", only_pos = TRUE)
# 
DoHeatmap(liver, features = unique(functions.markers$gene), assay = "function_n", angle = 90) + NoLegend()


liver_func.markers <- FindAllMarkers(liver, assay = "function_n", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
liver_func.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> res_cluster_top_gene

top_n_function_show <- liver_func.markers %>% group_by(cluster) %>% top_n(n = t_n, wt = avg_log2FC)
#write top n variable for each cluster


DoHeatmap(liver, features = unique(res_cluster_top_gene$gene), assay = "function_n", angle = 90) + NoLegend()

top_n_each_cluster <- paste(reactome_out_dir,"reactome_mm_top_n_variable_funciotn_for_each_cluster.csv",sep = "/")
write.csv(top_n_function_show, top_n_each_cluster, row.names = TRUE)


cluster_marker_gene <- unique(res_cluster_top_gene$gene)

for (ss in cluster_marker_gene){
  print(FeaturePlot(liver, features = ss, label = TRUE))
  tem <- readline("press any key to continue!")
}

for (ss in cluster_marker_gene){
  print(RidgePlot(liver, features = ss))
  tem <- readline("press any key to continue!")
}


#write top n variable for each cluster
cluster_all <- paste(reactome_out_dir, "reactome_mm_cluster_all.csv", sep = "/")
write.csv(liver@meta.data$seurat_clusters, cluster_all, row.names = FALSE)

# write cluster information for each clusters
cluster_marker_source <- paste(reactome_out_dir,"reactome_mm_cluster_marker_source.csv", sep = "/")
write.csv(liver.markers,cluster_marker_source)

# write cluster markers
cluster_marker_function <- paste(reactome_out_dir, "reactome_mm_cluster_marker_function.csv", sep = "/")
write.csv(liver_func.markers,cluster_marker_function)


#replot
cell_types_fn <- paste(data_dir,"cell_type.tsv",sep = "/")
cell_types <- read.table(cell_types_fn, stringsAsFactors = FALSE)
cell_types <- cell_types[,1]
cell_types_l <-  as.factor(cell_types)
dim(cell_types_l)
cell_types_l


liver@meta.data$cell_types <- cell_types_l
DimPlot(liver, reduction = "umap", group.by = "cell_types", label = TRUE) 

DimPlot(liver, reduction = "umap", group.by = "seurat_clusters", label = TRUE) 

top_n_show_ct <- liver.markers %>% group_by("cell_types") %>% top_n(n = t_n, wt = avg_log2FC)
DoHeatmap(liver, features = top_n_show_ct$gene) + NoLegend()




#############################################################################################
rm(list=ls())


# function on origin go mm!
#############################################################################################
# set directory
#############################################################################################
data_dir = "./data/human_liver"
out_dir = "./result/human_liver"

reactome_dir = "./data/gsea_cut_off/human_liver/human_liver_c2"
reactome_out_dir = "./result/gsea_cut_off/human_liver/human_liver_c2"

go_dir = "./data/gsea_cut_off/human_liver/human_liver_c5"
go_out_dir = "./result/gsea_cut_off/human_liver/human_liver_c5"

validation_dir = "./data/gsea_cut_off/human_liver/human_liver_c7"
validation_out_dir = "./result/gsea_cut_off/human_liver/human_liver_c7"
#############################################################################################


s_data_dir = data_dir

liver.data <- Read10X(data.dir = s_data_dir)

n_feature = 1
n_feature <- as.numeric(n_feature)

prject_name = "go"
liver <- CreateSeuratObject(counts = liver.data, project = "ORINGIN", min.cells = 3, min.features = n_feature)
liver
VlnPlot(liver, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)



nFeature_RNA_max <- 8000

nFeature_RNA_max <- as.numeric(nFeature_RNA_max)
scale_factor <- 35000
scale_factor <- as.numeric(scale_factor)


#normalizing the data
liver <- NormalizeData(liver, normalization.method = "LogNormalize", scale.factor = scale_factor)   
liver <- NormalizeData(liver)  

liver <- FindVariableFeatures(liver, selection.method = "vst", nfeatures = nFeature_RNA_max)    

t_n <- 5
top_n <- head(VariableFeatures(liver), t_n)


#scaling the data
all.genes <- rownames(liver)
liver <- ScaleData(liver, features = all.genes)

#run pca
liver <- RunPCA(liver, features = VariableFeatures(object = liver))    
print(liver[["pca"]], dims = 1:5, nfeatures = 5)


VizDimLoadings(liver, dims = 1:2, reduction = "pca")


DimPlot(liver, reduction = "pca")


ElbowPlot(liver)
umap_pca <- 15 
umap_pca <- as.numeric(umap_pca)

resolution_num = 0.5
resolution_num <- as.numeric(resolution_num)

liver <- FindNeighbors(liver, dims = 1:umap_pca)
liver <- FindClusters(liver, resolution = resolution_num)   

liver <- RunUMAP(liver, dims = 1:umap_pca)   

DimPlot(liver, reduction = "umap")    



liver.markers <- FindAllMarkers(liver, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
liver.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

DimPlot(liver, label = TRUE) + NoLegend()




#Add a function Seurat Object to plot!
liver_func.data <- Read10X(data.dir = go_dir)

n_feature_func = 2
n_feature_func <- as.numeric(n_feature_func)
func_obj_name = "function name!"
func_obj <- CreateSeuratObject(counts = liver_func.data, project = func_obj_name, min.cells = 3, min.features = n_feature_func)



VlnPlot(func_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)


liver[["function_n"]] = CreateAssayObject(liver_func.data)


func_feature_max <- 6000
func_feature_max <- as.numeric(func_feature_max)

func_scale_factor <- 20000
func_scale_factor <- as.numeric(func_scale_factor)

liver <- NormalizeData(liver, assay = "function_n", normalization.method = "LogNormalize",  scale.factor = scale_factor)

liver <- ScaleData(liver, assay = "function_n")

liver <- FindVariableFeatures(liver, assay = "function_n", selection.method = "vst", nfeatures = func_feature_max)

t_n_f <- 20
t_n_f <- as.numeric(t_n_f)
top_n_f <- head(VariableFeatures(liver, assay = "function_n"), t_n_f)


for (ss in top_n_f){
  print(FeaturePlot(liver, features = ss))
  tem <- readline("press any key to continue!")
}

for (ss in top_n_f){
  print(RidgePlot(liver, features = ss))
  tem <- readline("press any key to continue!")
}


# #Identify differentially expressed functions between clusters
# max_cells = 777
# liver.small <- subset(liver, downsample = max_cells)
# functions.markers <- FindAllMarkers(liver , assay = "function_n", only_pos = TRUE)
# DoHeatmap(liver.small, features = unique(functions.markers$gene), assay = "function_n", angle = 90) + NoLegend()


liver_func.markers <- FindAllMarkers(liver, assay = "function_n", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
liver_func.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> res_cluster_top_gene

top_n_function_show <- liver_func.markers %>% group_by(cluster) %>% top_n(n = t_n, wt = avg_log2FC)

DoHeatmap(liver, features = unique(res_cluster_top_gene$gene), assay = "function_n", angle = 90) + NoLegend()


#write top n variable for each cluster
top_n_each_cluster <- paste(go_out_dir,"go_mm_top_n_variable_funciotn_for_each_cluster.csv",sep = "/")
write.csv(top_n_function_show, top_n_each_cluster, row.names = TRUE)




cluster_marker_gene <- unique(res_cluster_top_gene$gene)

for (ss in cluster_marker_gene){
  print(FeaturePlot(liver, features = ss, label = TRUE))
  tem <- readline("press any key to continue!")
}

for (ss in cluster_marker_gene){
  print(RidgePlot(liver, features = ss))
  tem <- readline("press any key to continue!")
}



#write some information to output file!
cluster_all <- paste(go_out_dir, "go_mm_cluster_all.csv", sep = "/")

write.csv(liver@meta.data$seurat_clusters, cluster_all, row.names = FALSE)

# write cluster information for each clusters
cluster_marker_source <- paste(go_out_dir,"go_mm_cluster_marker_source.csv", sep = "/")

write.csv(liver.markers,cluster_marker_source)


#write cluster markers
cluster_marker_function <- paste(go_out_dir, "go_mm_cluster_marker_function.csv", sep = "/")
write.csv(liver_func.markers,cluster_marker_function)



#replot
cell_types_fn <- paste(data_dir,"cell_type.tsv",sep = "/")
cell_types <- read.table(cell_types_fn, stringsAsFactors = FALSE)
cell_types <- cell_types[,1]
cell_types_l <-  as.factor(cell_types)
dim(cell_types_l)
cell_types_l


liver@meta.data$cell_types <- cell_types_l
DimPlot(liver, reduction = "umap", group.by = "cell_types", label = TRUE) 

DimPlot(liver, reduction = "umap", group.by = "seurat_clusters", label = TRUE) 

top_n_show_ct <- liver_func.markers %>% group_by("cell_types") %>% top_n(n = t_n, wt = avg_log2FC)





rm(list=ls())
#############################################################################################
