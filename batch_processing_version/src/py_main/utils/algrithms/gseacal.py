import gseapy as gp


def calcu_gsea(gene_list, gmt, bg):
    gene_list = list(gene_list)
    enr2 = gp.enrichr(gene_list=gene_list,
                      # or gene_list=glist
                      description='test_name',
                      gene_sets=gmt.term_set,
                      background=bg,  # or the number of genes, e.g 20000
                      outdir=None,
                      cutoff=0.5,  # only used for testing.
                      verbose=True)
    return enr2.results

