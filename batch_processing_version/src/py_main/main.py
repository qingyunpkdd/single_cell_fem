from utils.algrithms.gseacal import calcu_gsea
from utils.data_processing.data_process import GetData, Gmt_stat
from utils.data_processing.out_put_file import Formater
from multiprocessing import JoinableQueue, Pool, Process, Manager
from copy import deepcopy
from collections import OrderedDict
from tqdm import tqdm
import time
from shutil import copyfile
import os
from utils.algrithms.hypergeom import calcu_hypergeom, hypergeom_single, hypergeom_row
import utils.algrithms.mat_multi_hypergeom as matmulhypergeom


def save_file(sample_gsea, formater, func_counts, sample_counts):
    formater.to_orig_format_matrix_new(sample_gsea)
    formater.to_orig_marmat_matrix_meta(func_counts, sample_counts)
    formater.to_orig_format_func()


def consummer(expr_dict, gsea_samples, lock):
    while True:

        try:
            lock.acquire()

        except:
            print("lock error!")
            break

        if len(expr_dict) == 0:
            print('terminate!')
            lock.release()
            break
        data_raw = deepcopy(expr_dict.popitem())
        data_name, data_x = data_raw

        data_data = data_x["data"]
        number_gene_set_s_rem = len(expr_dict) % data_x["numb_gene_set"]
        if number_gene_set_s_rem % 100 == 0:
            print(number_gene_set_s_rem)
        lock.release()

        try:
            tem_res = []
            expr_genes_d = list(data_data["expr_in_geneset"])

            for expr_in_geneset, expr_genes, colname_cell in zip(expr_genes_d,  data_data["expr_genes"], data_x["colname_cells"]):

                p_value = hypergeom_single(expr_in_geneset=expr_in_geneset,
                                           allgenes=data_data["number_all_genes"],
                                           _geneset=data_data["_gene_set"],
                                           expr_genes=expr_genes
                                           )
                tem_dic = OrderedDict()
                tem_dic["rowname_gmt"] = data_x["rowname_gmt"]
                tem_dic["colname_cell"] = colname_cell
                tem_dic["pvalue"] = p_value
                lock.acquire()
                gsea_samples[str(data_x["rowname_gmt"])+"_"+str(colname_cell)] = tem_dic
                lock.release()

        except:
            raise ("calcualte error!")
        finally:

            del tem_dic


def consummer_row(expr_dict, gsea_samples, lock):
    while True:

        try:
            lock.acquire()

        except:
            print("lock error!")
            break

        if len(expr_dict) == 0:
            print('terminate!')
            lock.release()
            break
        data_name, data_x = expr_dict.popitem()
        data_data = data_x["data"]
        number_cells_rem = len(expr_dict) % data_x["numb_gene_set"]
        if number_cells_rem % 100 == 0:
            print(number_cells_rem)
        lock.release()

        try:
            p_values = hypergeom_row(expr_in_genesets=data_data["expr_in_geneset"],
                                       allgenes=data_data["number_all_genes"],
                                       _geneset=data_data["_gene_set"],
                                       expr_genes=data_data["expr_genes"]
                                     )
            rowname_gmt = data_x["rowname_gmt"],
            colname_cell = data_x["colname_cells"]

            if p_values:
                lock.acquire()
                gsea_samples[data_name] = (rowname_gmt, colname_cell, p_values)
                lock.release()

        except:
            raise ("calcualte error!")



def producer(expr_dict, expr_obj, gmt_obj, lock):
    all_genes_list = matmulhypergeom.get_all_genes(expr_obj, gmt_obj)
    number_all_genes = len(all_genes_list)
    expmat, col_name, sample_expr_gene_num = matmulhypergeom.get_expmat(expr_obj, all_genes_list)
    gmtmat, row_name, pathway_gene_num = matmulhypergeom.get_gmtmat(gmt_obj, all_genes_list)
    intersect_matrix = matmulhypergeom.intersect_mat(expmat, gmtmat)
    # test_fewer = False
    # test_consumer = False
    for ip, rname in tqdm(enumerate(row_name)):
        _gene_set = pathway_gene_num[ip]
        expr_in_geneset_each_row = intersect_matrix[ip, ]
        data_row = {"rowname_gmt": rname,
                    "colname_cells": col_name,
                    "numb_cells": len(col_name),
                    "numb_gene_set": len(row_name),
                    "data": {"expr_in_geneset": list(expr_in_geneset_each_row),
                             "number_all_genes": number_all_genes,
                             "_gene_set": _gene_set,
                             "expr_genes": sample_expr_gene_num}}

        # if test_fewer and ip == 100:
        #     break

        expr_dict[ip] = data_row

    print("update finished")
    # consummer(expr_dict, gsea_samples, lock)
    # if test_consumer:
    #     consummer_row(expr_dict, gsea_samples, lock)


def pool_init(q, q_out):
    global queue  # make queue global in workers
    global queue_out
    queue = q
    queue_out = q_out


if __name__ == '__main__':

    batch_expr_dir = "../../data/expressions"
    batch_gmt_dir = "../../data/gene_sets"

    out_dir = "../../result/gsea_out"
    for expr_each in os.listdir(batch_expr_dir):
        if os.path.isdir(os.path.join(batch_expr_dir, expr_each)):
            print('----expr----' * 10)
            expr_dir = os.path.join(batch_expr_dir, expr_each)
            print(expr_dir)
            print('----expr----' * 10)

            for gmt_each in os.listdir(batch_gmt_dir):
                if os.path.isfile(os.path.join(batch_gmt_dir, gmt_each)):
                    print(os.path.join(batch_gmt_dir, gmt_each))

                    gmt_dir_file = os.path.join(batch_gmt_dir, gmt_each)
                    tem, *x = gmt_each.split('.')
                    each_out_dir = str(expr_each) + '_' + '' + str(tem)
                    each_out_dir = os.path.join(out_dir, each_out_dir)
                    os.mkdir(each_out_dir)
                    # expr_dir = "./data/hg19"
                    # gmt_dir_file = "./data/gmt/reactome.gmt"
                    expr_data = GetData(expr_dir)
                    gmt_data = Gmt_stat(gmt_dir_file)
                    expr_data.run_1()
                    gmt_data.get_gmtobj()
                    gmt_data.get_gmt_info()

                    sample_counts = len(expr_data.sample_geneset)
                    func_counts = len(gmt_data.term_names)

                    funclist = gmt_data.term_names
                    meta_data = expr_data.matrix_ex_meta_data_head

                    formater = Formater(meta_data, out_dir=each_out_dir, funclist=funclist)

                    mgr = Manager()
                    lock = mgr.Lock()
                    expr_dict = mgr.dict()

                    gsea_samples = mgr.dict()

                    n_consumer = 5
                    pool = Pool(n_consumer)

                    print('ready to run multiprocessing')
                    print('--' * 20)
                    print('uploading data...')
                    ti = time.time()
                    producer(expr_dict, expr_obj=expr_data, gmt_obj=gmt_data, lock=lock)
                    # prod = Process(target=producer, args=(expr_dict, expr_data, gmt_data,))
                    # prod.start()
                    # prod.join()
                    bg = len(expr_data.gene_list)
                    print('run calculator...')

                    for i in range(n_consumer):
                        pool.apply_async(consummer_row, args=(expr_dict, gsea_samples, lock,))
                    print("starting")
                    pool.close()
                    pool.join()

                    samples_gsea = matmulhypergeom.format_transfer_row(gsea_samples)

                    save_file(samples_gsea, formater, func_counts, sample_counts)
                    copyfile(os.path.join(expr_dir, 'barcodes.tsv'), os.path.join(each_out_dir, 'barcodes.tsv'))

                    pool.terminate()
                    # prod.terminate()
                    te = time.time()
                    print("speed: {}cells/(core*sec)".format(len(expr_data.sample_geneset)/((te-ti)*n_consumer)))
                    del lock
                    del expr_dict
                    del samples_gsea
                    del mgr
                    del expr_data
                    del gmt_data
                    del bg

                    del funclist
                    del sample_counts
                    del func_counts
                    del formater

                    del each_out_dir
                    del gmt_dir_file
                    del x
                    del meta_data
                    del tem

                    print("time used:", te - ti)
