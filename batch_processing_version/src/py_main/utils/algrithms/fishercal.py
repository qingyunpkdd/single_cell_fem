import numpy as np
import scipy as sp
import abc
import pandas as pd

from scipy import stats
from collections import OrderedDict, Counter

import time
import os
def calcu_fisher(gene_list, gmt, bg):
    ti = time.time()
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
            oddsratio, pvalue = stats.fisher_exact([[allgenes, _geneset], [expr_genes, expr_in_geneset]])
            #oddsratio means the expr_in_geneset is greater than theory(>1) or lettle than theory(<1)
            oddsratio_thred = 1
            pvalue_thred = 0.5
            if oddsratio >= oddsratio_thred and pvalue <= pvalue_thred:
                row = {'Term':s, 'Adjusted P-value':pvalue}
                res = res.append(row, ignore_index=True)
            #res.loc[index] = [term_name, pvalue]
    te = time.time()
    tm = "pid:" + str(os.getpid()) + "-|-time:" + str(te - ti)
    print(tm)
    return res

def calcu_fisherv2(gene_list, gmt_data_dict, bg):
    ti = time.time()
    min_count =1
    gene_list_n = len(gene_list)
    allgenes = bg
    _intersection = []

    _term_number = gmt_data_dict['term_number']
    res = pd.DataFrame(columns=['Term', 'Adjusted P-value'])

    oddsratio_thred = 1
    pvalue_thred = 0.5

    for gene in gene_list:
        if gene in gmt_data_dict['gene_term_dic']:
            tem = gmt_data_dict['gene_term_dic'][gene]
            if not _intersection:
                _intersection = tem
            else:
                _intersection.extend(tem)
    term_count = Counter(_intersection)
    for _term, _count in term_count.items():
        if _count >=min_count:
            term_n = _term_number[_term]
            term_genelist_n = term_count[_term]
            oddsratio, pvalue = stats.fisher_exact([[allgenes, term_n], [gene_list_n, term_genelist_n]])


            if oddsratio >= oddsratio_thred and pvalue <= pvalue_thred:
                row = {'Term': _term, 'Adjusted P-value': pvalue}
                res = res.append(row, ignore_index=True)
    te = time.time()
    tm ="pid:" + str(os.getpid()) + "-|-time:" + str(te-ti)
    print(tm)
    return res

