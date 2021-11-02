# How to use

## require

```python
pyqt5==5.15.1
scipy==1.5.2
numpy==1.19.1
pandas==1.1.1
tqdm==4.49.0
```



## GUI version

```python

# go to the gui__versioin\src\py_main directory
cd ./gui__versioin/src/py_main
python main.py


```



![image-20201202143317111](./image-20201202143317111.png)



expression data directory is the 10x format single cell expression data directory.

gmt data directory is the gene set data that similar to the misgdb format data directory.

output directory is the function expression data directory.





## Batch processing version

```python

# go to the \batch_processing_version\src\py_main
cd ./batch_processing_version/src/py_main
python main.py


# the default data and result directory are below:
batch_expr_dir = "../../data/expressions"
batch_gmt_dir = "../../data/gene_sets"
out_dir = "../../result/gsea_out"


```

## R and python script to reproduce three data sets (PBMC, liver, pancreas)

All data processing codes in Seurat are placed in the script\R_script directory.



## data download

Download link for all data:

PBMC dataset：
https://www.10xgenomics.com/resources/datasets
https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

reference :
10X Genomics official website (https://www.10xgenomics.com/resources/datasets)
Butler A, Hoffman P, Smibert P, Papalexi E, Satija R. 2018. Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nature Biotechnology. DOI: 10.1038/nbt.4096.



liver dataset：

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81252

reference :

Camp, J. G. et al. Multilineage communication regulates human liver bud development from pluripotency. Nature 546, 533–538 (2017)



pancreas dataset：

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85241

reference :

Muraro, M. J. et al. A Single-Cell Transcriptome Atlas of the Human Pancreas. Cell Syst 3, 385–394.e3 (2016)

