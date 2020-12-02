import pandas as pd

from scipy import stats
from collections import OrderedDict, Counter

import time
import os
def calcu_hypergeom(gene_list, gmt, bg):
    #ti = time.time()
    allgenes = bg
    _geneset = len(gmt.term_set)
    min_count = 1

    res = pd.DataFrame(columns=['Term','Adjusted P-value'])
    index = 0
    low_str = False
    subsets = sorted(gmt.term_set.keys())
    for s in subsets:
        # if low_str:
        #     geneset = [str.lower(gene) for gene in geneset]
        #     gene_list = [str.lower(gene) for gene in gene_list]
        category = gmt.term_set.get(s)
        _geneset = len(category)
        expr_genes = len(gene_list)
        hits = category.intersection(gene_list)
        expr_in_geneset = len(hits)
        if expr_in_geneset > min_count:
            #oddsratio, pvalue = stats.fisher_exact([[allgenes, _geneset], [expr_genes, expr_in_geneset]])

            pvalue = stats.hypergeom.sf(expr_in_geneset-1, allgenes, _geneset, expr_genes)
            #oddsratio means the expr_in_geneset is greater than theory(>1) or lettle than theory(<1)
            pvalue_thred = 0.5
            if pvalue <= pvalue_thred:
                row = {'Term':s, 'Adjusted P-value':pvalue}
                res = res.append(row, ignore_index=True)
            #res.loc[index] = [term_name, pvalue]
    #te = time.time()
    #tm = "pid:" + str(os.getpid()) + "-|-time:" + str(te - ti)
    #print(tm)
    return res

def hypergeom_single(expr_in_geneset, allgenes, _geneset, expr_genes, min_count=1):
    if expr_in_geneset <= min_count:
        pvalue = 1
    else:
        pvalue = stats.hypergeom.sf(expr_in_geneset-1, allgenes, _geneset, expr_genes)
    return pvalue

def hypergeom_row(expr_in_genesets, allgenes, _geneset, expr_genes, min_count=2):
    each_gmt_pvalues = []
    for ih, value in enumerate(expr_in_genesets):
        expr_in_geneset = value
        expr_gene = expr_genes[ih]
        if expr_in_geneset <= min_count:
            pvalue = 1
        else:
            pvalue = stats.hypergeom.sf(expr_in_geneset-1, allgenes, _geneset, expr_gene)
        each_gmt_pvalues.append(pvalue)
    return each_gmt_pvalues

