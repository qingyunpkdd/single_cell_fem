# set data and output directory
#############################################################################################
data_dir = "./data/hg19"
out_dir = "./result/hg19"

reactome_dir = "./data/gsea_cut_off/hg19/hg19_c2"
reactome_out_dir = "./result/gsea_cut_off/hg19/hg19_c2"

go_dir = "./data/gsea_cut_off/hg19/hg19_c5"
go_out_dir = "./result/gsea_cut_off/hg19/hg19_c5"

validation_dir = "./data/gsea_cut_off/hg19/hg19_c7"
validation_out_dir = "./result/gsea_cut_off/hg19/hg19_c7"
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
pbmc.data <- Read10X(data.dir = data_dir)
print(data_dir)
n_feature <- 1
prject_name = "origin"

pbmc <- CreateSeuratObject(counts = pbmc.data, project = prject_name, min.cells = 3, min.features = n_feature)
pbmc


nFeature_RNA_max <- 3000
nFeature_RNA_max <- as.numeric(nFeature_RNA_max)
scale_factor <- 10000
scale_factor <- as.numeric(scale_factor)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#normalizing the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = scale_factor)  
pbmc <- NormalizeData(pbmc)


pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = nFeature_RNA_max)    

t_n <- 5
t_n <- as.numeric(t_n)
top_n <- head(VariableFeatures(pbmc), t_n)



#scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#run pca
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))    
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5) 

umap_pca <- 10
umap_pca <- as.numeric(umap_pca)

resolution_num = 0.5
resolution_num <- as.numeric(resolution_num)

pbmc <- FindNeighbors(pbmc, dims = 1:umap_pca)
pbmc <- FindClusters(pbmc, resolution = resolution_num)    

pbmc <- RunUMAP(pbmc, dims = 1:umap_pca)    


# find all top n markers for each clusters
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)    
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


# plot top n markers for each clusters
for (ss in top_n){
  print(VlnPlot(pbmc, features = ss))
  tem <- readline("press any key to continue!")
}    
for (ss in top_n){
  print(FeaturePlot(pbmc, features = ss))
  tem <- readline("press any key to continue!")
}




#write top n variable for each cluster
top_n_show <- pbmc.markers %>% group_by(cluster) %>% top_n(n = t_n, wt = avg_log2FC)
top_n_each_cluster <- paste(out_dir,"top_n_variable_for_each_cluster.csv",sep = "/")
write.csv(top_n_show, top_n_each_cluster, row.names = TRUE)



# write cluster information for each clusters
cluster_all <- paste(out_dir,"cluster_all.csv",sep = "/")
write.csv(pbmc@meta.data$seurat_clusters, cluster_all, row.names = FALSE)

#write cluster markers
cluster_marker_all <- paste(out_dir,"cluster_marker.csv", sep = "/")
write.csv(pbmc.markers,cluster_marker_all)


# plot UMAP with cell type annotation
new.cluster.ids <- c("Naive_CD4_T", "CD14+_Mono", "Memory_CD4_T", "CD8_T", "B", "FCGR3A+_Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# write cell type information
cell_type_anno <- paste(out_dir,"cell_type.tsv", sep = "/")
write.csv(as.vector(pbmc@active.ident), cell_type_anno, row.names = FALSE, quote=F)

#############################################################################################



rm(list=ls())
# function run reactome
#############################################################################################

# set data and output directory
#############################################################################################
data_dir = "./data/hg19"
out_dir = "./result/hg19"

reactome_dir = "./data/gsea_cut_off/hg19/hg19_c2"
reactome_out_dir = "./result/gsea_cut_off/hg19/hg19_c2"

go_dir = "./data/gsea_cut_off/hg19/hg19_c5"
go_out_dir = "./result/gsea_cut_off/hg19/hg19_c5"

validation_dir = "./data/gsea_cut_off/hg19/hg19_c7"
validation_out_dir = "./result/gsea_cut_off/hg19/hg19_c7"
#############################################################################################





pbmc.data <- Read10X(data.dir = reactome_dir)
print(reactome_dir)
n_feature = 1
n_feature <- as.numeric(n_feature)
prject_name = "reactome"

pbmc <- CreateSeuratObject(counts = pbmc.data, project = prject_name, min.cells = 3, min.features = n_feature)
pbmc

nFeature_RNA_max <- 1400
nFeature_RNA_max <- as.numeric(nFeature_RNA_max)
scale_factor <- 8000
scale_factor <- as.numeric(scale_factor)


#normalizing the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = scale_factor)    
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = nFeature_RNA_max)    

t_n <- 5
t_n <- as.numeric(t_n)
top_n <- head(VariableFeatures(pbmc), t_n)


#scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


#run pca
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))    
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5) 

umap_pca <- 20
umap_pca <- as.numeric(umap_pca)


resolution_num = 0.6
resolution_num <- as.numeric(resolution_num)

pbmc <- FindNeighbors(pbmc, dims = 1:umap_pca)
pbmc <- FindClusters(pbmc, resolution = resolution_num)    
pbmc <- RunUMAP(pbmc, dims = 1:umap_pca)    

DimPlot(pbmc, reduction = "umap")    

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
# 
# for (ss in top_n){
#   print(VlnPlot(pbmc, features = ss))
#   tem <- readline("press any key to continue!")
# }    
# for (ss in top_n){
#   print(FeaturePlot(pbmc, features = ss))
#   tem <- readline("press any key to continue!")
# }

top_n_show <- pbmc.markers %>% group_by(cluster) %>% top_n(n = t_n, wt = avg_log2FC)
DoHeatmap(pbmc, features = top_n_show$gene) + NoLegend()


top_n_each_cluster_show <- unique(top_n_show$gene)
for (ss in top_n_each_cluster_show){
  print(VlnPlot(pbmc, features = ss))
  tem <- readline("press any key to continue!")
}

for (ss in top_n_each_cluster_show){
  print(FeaturePlot(pbmc, features = ss))
  tem <- readline("press any key to continue!")
}



#write top n variable for each cluster
top_n_each_cluster <- paste(reactome_out_dir,"reactome_top_n_variable_for_each_cluster.csv",sep = "/")
write.csv(top_n_show, top_n_each_cluster, row.names = TRUE)


# write cluster information for each clusters
cluster_all <- paste(reactome_out_dir,"reactome_cluster_all.csv",sep = "/")
write.csv(pbmc@meta.data$seurat_clusters, cluster_all, row.names = FALSE)

#write cluster markers
cluster_marker_all <- paste(reactome_out_dir,"reactome_cluster_marker.csv", sep = "/")
write.csv(pbmc.markers,cluster_marker_all)



#replot
cell_types_fn <- paste(data_dir,"cell_type.tsv",sep = "/")
cell_types <- read.table(cell_types_fn, stringsAsFactors = FALSE)
cell_types <- cell_types[,1]
cell_types_l <-  as.factor(cell_types)
dim(cell_types_l)
cell_types_l

# plot all umap
pbmc@meta.data$cell_types <- cell_types_l
DimPlot(pbmc, reduction = "umap", group.by = "cell_types", label = TRUE)
DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", label = TRUE)




#heatmap
library(ggplot2)
data_heatmap <- pbmc@assays$RNA@scale.data
data_row_name <- pbmc@assays$RNA@data@Dimnames[[1]]
data_col_name <- pbmc@assays$RNA@data@Dimnames[[2]]
row.names(data_heatmap) <- data_row_name
colnames(data_heatmap) <- data_col_name

heatmap(data_heatmap,cexRow=0.5)


#############################################################################################

rm(list=ls())





# function run go
#############################################################################################

# set data and output directory
#############################################################################################
data_dir = "./data/hg19"
out_dir = "./result/hg19"

reactome_dir = "./data/gsea_cut_off/hg19/hg19_c2"
reactome_out_dir = "./result/gsea_cut_off/hg19/hg19_c2"

go_dir = "./data/gsea_cut_off/hg19/hg19_c5"
go_out_dir = "./result/gsea_cut_off/hg19/hg19_c5"

validation_dir = "./data/gsea_cut_off/hg19/hg19_c7"
validation_out_dir = "./result/gsea_cut_off/hg19/hg19_c7"
#############################################################################################



pbmc.data <- Read10X(data.dir = go_dir)
print(go_dir)
n_feature = 1
n_feature <- as.numeric(n_feature)
prject_name = "go"
pbmc <- CreateSeuratObject(counts = pbmc.data, project = prject_name, min.cells = 3, min.features = n_feature)
pbmc



nFeature_RNA_max <- 7600
nFeature_RNA_max <- as.numeric(nFeature_RNA_max)
scale_factor <- 30000
scale_factor <- as.numeric(scale_factor)



#normalizing the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = scale_factor)    
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = nFeature_RNA_max)    

t_n <- 5
t_n <- as.numeric(t_n)
top_n <- head(VariableFeatures(pbmc), t_n)


#scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#run pca
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))    
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5) 


umap_pca <- 15
umap_pca <- as.numeric(umap_pca)


resolution_num = 0.8
resolution_num <- as.numeric(resolution_num)

pbmc <- FindNeighbors(pbmc, dims = 1:umap_pca)
pbmc <- FindClusters(pbmc, resolution = resolution_num)    

pbmc <- RunUMAP(pbmc, dims = 1:umap_pca)    


# DimPlot(pbmc, reduction = "umap")    
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


for (ss in top_n){
  print(VlnPlot(pbmc, features = ss))
  tem <- readline("press any key to continue!")
}    
for (ss in top_n){
  print(FeaturePlot(pbmc, features = ss))
  tem <- readline("press any key to continue!")
  
}


# plot each unique genes!
top_n_show <- pbmc.markers %>% group_by(cluster) %>% top_n(n = t_n, wt = avg_log2FC)
# DoHeatmap(pbmc, features = top_n_show$gene) + NoLegend()


top_n_each_cluster_show <- unique(top_n_show$gene)
for (ss in top_n_each_cluster_show){
  print(VlnPlot(pbmc, features = ss))
  tem <- readline("press any key to continue!")
}  

for (ss in top_n_each_cluster_show){
  print(FeaturePlot(pbmc, features = ss))
  tem <- readline("press any key to continue!")
  
}



#write top n variable for each cluster
top_n_each_cluster <- paste(go_out_dir,"go_top_n_variable_for_each_cluster.csv",sep = "/")
write.csv(top_n_show, top_n_each_cluster, row.names = TRUE)



#write some information to output file!
cluster_all <- paste(go_out_dir,"go_cluster_all.csv",sep = "/")
write.csv(pbmc@meta.data$seurat_clusters, cluster_all, row.names = FALSE)


cluster_marker_all <- paste(go_out_dir,"go_cluster_marker.csv", sep = "/")
write.csv(pbmc.markers,cluster_marker_all)


#replot
cell_types_fn <- paste(data_dir,"cell_type.tsv",sep = "/")
cell_types <- read.table(cell_types_fn, stringsAsFactors = FALSE)
cell_types <- cell_types[,1]
cell_types_l <-  as.factor(cell_types)
length(cell_types_l)
cell_types_l


pbmc@meta.data$cell_types <- cell_types_l
DimPlot(pbmc, reduction = "umap", group.by = "cell_types", label = TRUE) 

DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

top_n_show_ct <- pbmc.markers %>% group_by("cell_types") %>% top_n(n = t_n, wt = avg_log2FC)
DoHeatmap(pbmc, features = top_n_show_ct$gene) + NoLegend()



#heatmap
library(ggplot2)
data_heatmap <- pbmc@assays$RNA@scale.data
data_row_name <- pbmc@assays$RNA@data@Dimnames[[1]]
data_col_name <- pbmc@assays$RNA@data@Dimnames[[2]]
row.names(data_heatmap) <- data_row_name
colnames(data_heatmap) <- data_col_name

heatmap(data_heatmap,cexRow=0.5)

#############################################################################################


rm(list=ls())

# function on origin Reactome mm !
#############################################################################################

# set data and output directory
#############################################################################################
data_dir = "./data/hg19"
out_dir = "./result/hg19"

reactome_dir = "./data/gsea_cut_off/hg19/hg19_c2"
reactome_out_dir = "./result/gsea_cut_off/hg19/hg19_c2"

go_dir = "./data/gsea_cut_off/hg19/hg19_c5"
go_out_dir = "./result/gsea_cut_off/hg19/hg19_c5"

validation_dir = "./data/gsea_cut_off/hg19/hg19_c7"
validation_out_dir = "./result/gsea_cut_off/hg19/hg19_c7"
#############################################################################################

s_data_dir <- data_dir

pbmc.data <- Read10X(data.dir = s_data_dir)

n_feature = 1
n_feature <- as.numeric(n_feature)

prject_name = "REACTOME"
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "ORINGIN", min.cells = 3, min.features = n_feature)

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)


nFeature_RNA_max <- 3000
nFeature_RNA_max <- as.numeric(nFeature_RNA_max)
scale_factor <- 10000
scale_factor <- as.numeric(scale_factor)


#normalizing the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = scale_factor)
pbmc <- NormalizeData(pbmc)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = nFeature_RNA_max)    

t_n <- 5
t_n <- as.numeric(t_n)
top_n <- head(VariableFeatures(pbmc), t_n)

plot2_1 <- VariableFeaturePlot(pbmc)
plot2_2 <- LabelPoints(plot = plot2_1, points = top_n, repel = TRUE)
CombinePlots(plots = list(plot2_1, plot2_2))     



#scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#run pca
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))    
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")


DimPlot(pbmc, reduction = "pca")


ElbowPlot(pbmc)
umap_pca <- 10
umap_pca <- as.numeric(umap_pca)

resolution_num = 0.5
resolution_num <- as.numeric(resolution_num)

pbmc <- FindNeighbors(pbmc, dims = 1:umap_pca)
pbmc <- FindClusters(pbmc, resolution = resolution_num)   

pbmc <- RunUMAP(pbmc, dims = 1:umap_pca)   

DimPlot(pbmc, reduction = "umap")    



pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

DimPlot(pbmc, label = TRUE) + NoLegend()




#Add a function Seurat Object to plot!
pbmc_func.data <- Read10X(data.dir = reactome_dir)

n_feature_func = 2
n_feature_func <- as.numeric(n_feature_func)
func_obj_name = "function name!"
func_obj <- CreateSeuratObject(counts = pbmc_func.data, project = func_obj_name, min.cells = 3, min.features = n_feature_func)


pbmc[["function_n"]] = CreateAssayObject(pbmc_func.data)


func_feature_max <- 1400

func_feature_max <- as.numeric(func_feature_max)

func_scale_factor <- 6000
func_scale_factor <- as.numeric(func_scale_factor)

pbmc <- NormalizeData(pbmc, assay = "function_n", normalization.method = "LogNormalize",  scale.factor = scale_factor)

pbmc <- ScaleData(pbmc, assay = "function_n")

pbmc <- FindVariableFeatures(pbmc, assay = "function_n", selection.method = "vst", nfeatures = func_feature_max)

t_n_f <- 20
t_n_f <- as.numeric(t_n_f)
top_n_f <- head(VariableFeatures(pbmc, assay = "function_n"), t_n_f)

special_function <- "reactome-response-to-elevated-platelet-cytosolic-ca2plus"
print(FeaturePlot(pbmc, features = special_function))
print(RidgePlot(pbmc, features = special_function))

for (ss in top_n_f){
  print(FeaturePlot(pbmc, features = ss))
  tem <- readline("press any key to continue!")
}

for (ss in top_n_f){
  print(RidgePlot(pbmc, features = ss))
  tem <- readline("press any key to continue!")
}


#Identify differentially expressed functions between clusters

max_cells = 777

pbmc.small <- subset(pbmc, downsample = max_cells)

functions.markers <- FindAllMarkers(pbmc , assay = "function_n", only_pos = TRUE)

DoHeatmap(pbmc, features = unique(functions.markers$gene), assay = "function_n", angle = 90) + NoLegend()


pbmc_func.markers <- FindAllMarkers(pbmc, assay = "function_n", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc_func.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> res_cluster_top_gene

top_n_function_show <- pbmc_func.markers %>% group_by(cluster) %>% top_n(n = t_n, wt = avg_log2FC)
#write top n variable for each cluster


DoHeatmap(pbmc, features = unique(res_cluster_top_gene$gene), assay = "function_n", angle = 90) + NoLegend()

top_n_each_cluster <- paste(reactome_out_dir,"reactome_mm_top_n_variable_funciotn_for_each_cluster.csv",sep = "/")
write.csv(top_n_function_show, top_n_each_cluster, row.names = TRUE)


cluster_marker_gene <- unique(res_cluster_top_gene$gene)

for (ss in cluster_marker_gene){
  print(FeaturePlot(pbmc, features = ss, label = TRUE))
  tem <- readline("press any key to continue!")
}

for (ss in cluster_marker_gene){
  print(RidgePlot(pbmc, features = ss))
  tem <- readline("press any key to continue!")
}


# write cluster information for each clusters
cluster_all <- paste(reactome_out_dir, "reactome_mm_cluster_all.csv", sep = "/")
write.csv(pbmc@meta.data$seurat_clusters, cluster_all, row.names = FALSE)

#write cluster markers(genes)
cluster_marker_source <- paste(reactome_out_dir,"reactome_mm_cluster_marker_source.csv", sep = "/")
write.csv(pbmc.markers,cluster_marker_source)

#write cluster markers(function)
cluster_marker_function <- paste(reactome_out_dir, "reactome_mm_cluster_marker_function.csv", sep = "/")
write.csv(pbmc_func.markers,cluster_marker_function)


#replot
cell_types_fn <- paste(data_dir,"cell_type.tsv",sep = "/")
cell_types <- read.table(cell_types_fn, stringsAsFactors = FALSE)
cell_types <- cell_types[,1]
cell_types_l <-  as.factor(cell_types)
length(cell_types_l)
cell_types_l


pbmc@meta.data$cell_types <- cell_types_l
DimPlot(pbmc, reduction = "umap", group.by = "cell_types", label = TRUE) 

DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", label = TRUE) 



#############################################################################################
rm(list=ls())


# function on origin go mm!
#############################################################################################

# set data and output directory
#############################################################################################
data_dir = "./data/hg19"
out_dir = "./result/hg19"

reactome_dir = "./data/gsea_cut_off/hg19/hg19_c2"
reactome_out_dir = "./result/gsea_cut_off/hg19/hg19_c2"

go_dir = "./data/gsea_cut_off/hg19/hg19_c5"
go_out_dir = "./result/gsea_cut_off/hg19/hg19_c5"

validation_dir = "./data/gsea_cut_off/hg19/hg19_c7"
validation_out_dir = "./result/gsea_cut_off/hg19/hg19_c7"
#############################################################################################



s_data_dir <- data_dir


pbmc.data <- Read10X(data.dir = s_data_dir)

n_feature = 1
n_feature <- as.numeric(n_feature)

prject_name = "go"
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "ORINGIN", min.cells = 3, min.features = n_feature)
pbmc

nFeature_RNA_max <- 3000

nFeature_RNA_max <- as.numeric(nFeature_RNA_max)
scale_factor <- 10000
scale_factor <- as.numeric(scale_factor)


#normalizing the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = scale_factor)   
pbmc <- NormalizeData(pbmc)  

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = nFeature_RNA_max)    

t_n <- 5
t_n <- as.numeric(t_n)
top_n <- head(VariableFeatures(pbmc), t_n)


#scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#run pca
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))    
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)


# umap_pca <- readline("please input the first n dim to run the UMAP/T-SNE algrithm! :") 
umap_pca <- 10
umap_pca <- as.numeric(umap_pca)

# resolution_num = readline("please input the cluster factor to run cluster (ref:0.4-1.2,bigger number representate the greater number of cluster) :")
resolution_num = 0.5
resolution_num <- as.numeric(resolution_num)

pbmc <- FindNeighbors(pbmc, dims = 1:umap_pca)
pbmc <- FindClusters(pbmc, resolution = resolution_num)   

pbmc <- RunUMAP(pbmc, dims = 1:umap_pca)   

DimPlot(pbmc, reduction = "umap")

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

DimPlot(pbmc, label = TRUE) + NoLegend()


#Add a function Seurat Object to plot!
pbmc_func.data <- Read10X(data.dir = go_dir)

n_feature_func = 2
n_feature_func <- as.numeric(n_feature_func)
func_obj_name = "function name!"
func_obj <- CreateSeuratObject(counts = pbmc_func.data, project = func_obj_name, min.cells = 3, min.features = n_feature_func)



VlnPlot(func_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)


pbmc[["function_n"]] = CreateAssayObject(pbmc_func.data)


func_feature_max <- 6500

func_feature_max <- as.numeric(func_feature_max)

func_scale_factor <- 30000
func_scale_factor <- as.numeric(func_scale_factor)

pbmc <- NormalizeData(pbmc, assay = "function_n", normalization.method = "LogNormalize",  scale.factor = scale_factor)

pbmc <- ScaleData(pbmc, assay = "function_n")

pbmc <- FindVariableFeatures(pbmc, assay = "function_n", selection.method = "vst", nfeatures = func_feature_max)

t_n_f <- 20
t_n_f <- as.numeric(t_n_f)
top_n_f <- head(VariableFeatures(pbmc, assay = "function_n"), t_n_f)


for (ss in top_n_f){
  print(FeaturePlot(pbmc, features = ss))
  tem <- readline("press any key to continue!")
}

for (ss in top_n_f){
  print(RidgePlot(pbmc, features = ss))
  tem <- readline("press any key to continue!")
}


#Identify differentially expressed functions between clusters
# 
# max_cells = 2700
# 
# pbmc.small <- subset(pbmc, downsample = max_cells)
# 
# functions.markers <- FindAllMarkers(pbmc , assay = "function_n", only_pos = TRUE)
# 
# DoHeatmap(pbmc.small, features = unique(functions.markers$gene), assay = "function_n", angle = 90) + NoLegend()



pbmc_func.markers <- FindAllMarkers(pbmc, assay = "function_n", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc_func.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> res_cluster_top_gene

top_n_function_show <- pbmc_func.markers %>% group_by(cluster) %>% top_n(n = t_n, wt = avg_log2FC)

DoHeatmap(pbmc, features = unique(res_cluster_top_gene$gene), assay = "function_n", angle = 90) + NoLegend()


#write top n variable for each cluster



top_n_each_cluster <- paste(go_out_dir,"go_mm_top_n_variable_funciotn_for_each_cluster.csv",sep = "/")
write.csv(top_n_function_show, top_n_each_cluster, row.names = TRUE)




cluster_marker_gene <- unique(res_cluster_top_gene$gene)

for (ss in cluster_marker_gene){
  print(FeaturePlot(pbmc, features = ss, label = TRUE))
  tem <- readline("press any key to continue!")
}

for (ss in cluster_marker_gene){
  print(RidgePlot(pbmc, features = ss))
  tem <- readline("press any key to continue!")
}



# write cluster information for each clusters
cluster_all <- paste(go_out_dir, "go_mm_cluster_all.csv", sep = "/")
write.csv(pbmc@meta.data$seurat_clusters, cluster_all, row.names = FALSE)

#write cluster markers(genes)
cluster_marker_source <- paste(go_out_dir,"go_mm_cluster_marker_source.csv", sep = "/")
write.csv(pbmc.markers,cluster_marker_source)

#write cluster markers(function)
cluster_marker_function <- paste(go_out_dir, "go_mm_cluster_marker_function.csv", sep = "/")
write.csv(pbmc_func.markers,cluster_marker_function)



#replot
cell_types_fn <- paste(data_dir,"cell_type.tsv",sep = "/")
cell_types <- read.table(cell_types_fn, stringsAsFactors = FALSE)
cell_types <- cell_types[,1]
cell_types_l <-  as.factor(cell_types)
length(cell_types_l)
cell_types_l


pbmc@meta.data$cell_types <- cell_types_l
DimPlot(pbmc, reduction = "umap", group.by = "cell_types", label = TRUE) 
DimPlot(pbmc, reduction = "umap", group.by = "cell_types", label = FALSE) 

DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", label = TRUE) 
DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", label = FALSE) 

pbmc@meta.data$seurat_clusters <- cell_types_l
pbmc_func.markers <- FindAllMarkers(pbmc, assay = "function_n", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top_n_show_ct <- pbmc_func.markers %>% group_by("cell_types") %>% top_n(n = t_n, wt = avg_log2FC)
DoHeatmap(pbmc, features = top_n_show_ct$gene) + NoLegend()






rm(list=ls())
#############################################################################################










# set data and output directory
#############################################################################################
data_dir = "./data/hg19"
out_dir = "./result/hg19"

reactome_dir = "./data/gsea_cut_off/hg19/hg19_c2"
reactome_out_dir = "./result/gsea_cut_off/hg19/hg19_c2"

go_dir = "./data/gsea_cut_off/hg19/hg19_c5"
go_out_dir = "./result/gsea_cut_off/hg19/hg19_c5"

validation_dir = "./data/gsea_cut_off/hg19/hg19_c7"
validation_out_dir = "./result/gsea_cut_off/hg19/hg19_c7"
#############################################################################################


validation
#############################################################################################
s_data_dir <- data_dir


pbmc.data <- Read10X(data.dir = s_data_dir)

n_feature = 1
n_feature <- as.numeric(n_feature)

prject_name = "go"
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "ORINGIN", min.cells = 3, min.features = n_feature)
pbmc

nFeature_RNA_max <- 3000

nFeature_RNA_max <- as.numeric(nFeature_RNA_max)
scale_factor <- 10000
scale_factor <- as.numeric(scale_factor)


#normalizing the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = scale_factor)   
pbmc <- NormalizeData(pbmc)  

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = nFeature_RNA_max)    

t_n <- 5
t_n <- as.numeric(t_n)
top_n <- head(VariableFeatures(pbmc), t_n)


#scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#run pca
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))    
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)


# umap_pca <- readline("please input the first n dim to run the UMAP/T-SNE algrithm! :") 
umap_pca <- 10
umap_pca <- as.numeric(umap_pca)

# resolution_num = readline("please input the cluster factor to run cluster (ref:0.4-1.2,bigger number representate the greater number of cluster) :")
resolution_num = 0.5
resolution_num <- as.numeric(resolution_num)

pbmc <- FindNeighbors(pbmc, dims = 1:umap_pca)
pbmc <- FindClusters(pbmc, resolution = resolution_num)   

pbmc <- RunUMAP(pbmc, dims = 1:umap_pca)   

DimPlot(pbmc, reduction = "umap")

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

DimPlot(pbmc, label = TRUE) + NoLegend()


#Add a function Seurat Object to plot!
pbmc_func.data <- Read10X(data.dir = validation_dir)

n_feature_func = 2
n_feature_func <- as.numeric(n_feature_func)
func_obj_name = "function name!"
func_obj <- CreateSeuratObject(counts = pbmc_func.data, project = func_obj_name, min.cells = 3, min.features = n_feature_func)



VlnPlot(func_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)


pbmc[["function_n"]] = CreateAssayObject(pbmc_func.data)


func_feature_max <- 6500

func_feature_max <- as.numeric(func_feature_max)

func_scale_factor <- 30000
func_scale_factor <- as.numeric(func_scale_factor)

pbmc <- NormalizeData(pbmc, assay = "function_n", normalization.method = "LogNormalize",  scale.factor = scale_factor)

pbmc <- ScaleData(pbmc, assay = "function_n")

pbmc <- FindVariableFeatures(pbmc, assay = "function_n", selection.method = "vst", nfeatures = func_feature_max)

t_n_f <- 20
t_n_f <- as.numeric(t_n_f)
top_n_f <- head(VariableFeatures(pbmc, assay = "function_n"), t_n_f)


for (ss in top_n_f){
  print(FeaturePlot(pbmc, features = ss))
  tem <- readline("press any key to continue!")
}

for (ss in top_n_f){
  print(RidgePlot(pbmc, features = ss))
  tem <- readline("press any key to continue!")
}


#Identify differentially expressed functions between clusters

max_cells = 2700

pbmc.small <- subset(pbmc, downsample = max_cells)

functions.markers <- FindAllMarkers(pbmc , assay = "function_n", only_pos = TRUE)

DoHeatmap(pbmc.small, features = unique(functions.markers$gene), assay = "function_n", angle = 90) + NoLegend()

tem <- readline("press any key to continue!") 

pbmc_func.markers <- FindAllMarkers(pbmc, assay = "function_n", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc_func.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> res_cluster_top_gene

top_n_function_show <- pbmc_func.markers %>% group_by(cluster) %>% top_n(n = t_n, wt = avg_log2FC)

DoHeatmap(pbmc, features = unique(res_cluster_top_gene$gene), assay = "function_n", angle = 90) + NoLegend()


#write top n variable for each cluster

top_n_each_cluster <- paste(go_out_dir,"go_mm_top_n_variable_funciotn_for_each_cluster.csv",sep = "/")
write.csv(top_n_function_show, top_n_each_cluster, row.names = TRUE)

plot_markers <- c("gse10325-cd4-tcell-vs-myeloid-dn",
                  "gse11057-pbmc-vs-mem-cd4-tcell-up",
                  "gse22886-naive-bcell-vs-monocyte-up",
                  "gse22886-naive-tcell-vs-monocyte-up",
                  "gse26495-naive-vs-pd1low-cd8-tcell-up",
                  "gse7764-nkcell-vs-splenocyte-up"
)

eps("bac.txt")


cluster_marker_gene <- unique(res_cluster_top_gene$gene)

for (sc in plot_markers){
  print(FeaturePlot(pbmc, features = sc))
  plot_out_dir <- paste(validation_out_dir, paste(sc,".eps", sep = ""), sep = "/" )
  postscript(plot_out_dir)
  print(FeaturePlot(pbmc, features = sc))
  # postscript(plot_out_dir)
  dev.off()
  tem <- readline("press any key to continue!")
}

for (ss in cluster_marker_gene){
  print(RidgePlot(pbmc, features = ss))
  tem <- readline("press any key to continue!")
}



#write some information to output file!
cluster_all <- paste(validation_out_dir, "go_mm_cluster_all.csv", sep = "/")
write.csv(pbmc@meta.data$seurat_clusters, cluster_all, row.names = FALSE)

#write cluster markers(genes)
cluster_marker_source <- paste(validation_out_dir,"go_mm_cluster_marker_source.csv", sep = "/")
write.csv(pbmc.markers,cluster_marker_source)

#write cluster markers(function)
cluster_marker_function <- paste(validation_out_dir, "go_mm_cluster_marker_function.csv", sep = "/")
write.csv(pbmc_func.markers,cluster_marker_function)



#replot
cell_types_fn <- paste(data_dir,"cell_type.tsv",sep = "/")
cell_types <- read.table(cell_types_fn, stringsAsFactors = FALSE)
cell_types <- cell_types[,1]
cell_types_l <-  as.factor(cell_types)
dim(cell_types_l)
cell_types_l


pbmc@meta.data$cell_types <- cell_types_l
DimPlot(pbmc, reduction = "umap", group.by = "cell_types", label = TRUE) 
DimPlot(pbmc, reduction = "umap", group.by = "cell_types", label = FALSE) 

DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", label = TRUE) 
DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", label = FALSE) 



top_n_show_ct <- pbmc_func.markers %>% group_by("cell_types") %>% top_n(n = t_n, wt = avg_log2FC)
DoHeatmap(pbmc, features = top_n_show_ct$gene) + NoLegend()

