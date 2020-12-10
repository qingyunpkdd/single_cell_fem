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



![image-20201202143317111](。/image-20201202143317111.png)



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









