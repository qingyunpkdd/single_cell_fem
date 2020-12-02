import numpy as np
import pandas as pd
from .hypergeom import hypergeom_single
from collections import OrderedDict
from tqdm import tqdm

def get_all_genes(expr_obj, gmt_obj):
    expr_gene = expr_obj.gene_list
    gmt_gene = gmt_obj.term_all_gene
    all_genes = set(expr_gene)
    all_genes.update(set(gmt_gene))
    all_genes_res = []
    for gene in all_genes:
        all_genes_res.append(gene.lower())
    return list(set(all_genes_res))


def get_expmat(expr_obj, all_genes):
    num_row = len(all_genes)
    num_col = len(expr_obj.sample_geneset)
    expmat = np.zeros((num_row, num_col))
    all_gene_dic = {gene: i for i, gene in enumerate(all_genes)}
    for sid, geneset in expr_obj.sample_geneset.items():
        sid = int(sid)
        gene_set_each = [gene.lower() for gene in geneset]
        row_index = [all_gene_dic[gene] for gene in gene_set_each]
        # row_index = [all_genes.index(gene) for gene in gene_set_each]
        for r_idx in row_index:
            expmat[r_idx, sid] = 1

    sample_expr_gene_num = [len(geneset) for geneset in expr_obj.sample_geneset.values()]
    col_name = [sid for sid in expr_obj.sample_geneset.keys()]
    return expmat, col_name, sample_expr_gene_num


def get_gmtmat(gmt_obj, all_genes):
    num_col = len(all_genes)
    num_row = gmt_obj.term_quantity
    gmtmat = np.zeros((num_row, num_col))
    all_gene_dic = {gene: i for i, gene in enumerate(all_genes)}
    #row_name = [gmt_id for gmt_id in gmt_obj.term_set.keys()]
    row_name = []
    gmt_id = 0
    for gmt_name, geneset_i in gmt_obj.term_set.items():
        gene_set_each = [gene.lower() for gene in geneset_i]
        col_index = [all_gene_dic[gene] for gene in gene_set_each]
        row_name.append(gmt_name)
        for col_idx in col_index:
            gmtmat[gmt_id, col_idx] = 1
        gmt_id += 1
    pathway_gene_num = [len(geneset) for geneset in gmt_obj.term_set.values()]

    return gmtmat, row_name, pathway_gene_num


# exprmat: number of genes * number of cells
# gmtmat: number of pythway * number of genes
# res: number of pathway * numner of cells
def intersect_mat(exprmat, gmtmat):
    res = np.dot(gmtmat, exprmat)
    return res


def format_transfer(gsea_samples):
    samples_gsea = OrderedDict()
    samples_i = OrderedDict()
    res = pd.DataFrame(columns=['Term', 'Adjusted P-value'])
    for data_name, cont in gsea_samples.items():
        cell_id = cont["colname_cell"]
        gmt_name = cont["rowname_gmt"]
        if cell_id in samples_i:
            samples_i[cell_id].append((gmt_name, cont["pvalue"]))
        else:
            samples_i[cell_id] = [(gmt_name, cont["pvalue"])]

    # for cell_id, data_list in tqdm(samples_i.items()):
    #     res = pd.DataFrame(columns=['Term', 'Adjusted P-value'])
    #     for gmt_name, p_value in data_list:
    #         row = {'Term': gmt_name, 'Adjusted P-value': p_value}
    #         res = res.append(row, ignore_index=True)
    #     samples_gsea[cell_id] = res

    return samples_i


def format_transfer_row(gsea_samples):
    samples_i = OrderedDict()
    for data_name, cont in gsea_samples.items():
        rown = cont[0][0]
        for coln, pv in zip(cont[1], cont[2]):
            if coln in samples_i:
                samples_i[coln].append((rown, pv))
            else:
                samples_i[coln] = [(rown, pv)]

    return samples_i
