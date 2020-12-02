from utils.data_processing.data_process import GetData, Gmt_stat
from utils.data_processing.out_put_file import Formater
from multiprocessing import JoinableQueue, Pool, Process, Manager, freeze_support
from copy import deepcopy
from collections import OrderedDict
from tqdm import tqdm
import time
from shutil import copyfile
import os
from utils.algrithms.hypergeom import calcu_hypergeom, hypergeom_single, hypergeom_row
import utils.algrithms.mat_multi_hypergeom as matmulhypergeom
from utils.gui import main_win_setup
from PyQt5.QtWidgets import QApplication
import sys


def save_file(sample_gsea, formater, func_counts, sample_counts):
    formater.to_orig_format_matrix_new(sample_gsea)
    formater.to_orig_marmat_matrix_meta(func_counts, sample_counts)
    formater.to_orig_format_func()

def consummer_row(expr_dict, gsea_samples, lock, status_que=None):
    while True:

        try:
            lock.acquire()

        except:
            # print("lock error!")
            status_que.put({"print_out": "lock error!"})
            break

        if len(expr_dict) == 0:
            # print('terminate!')
            status_que.put({"print_out": "terminate"})
            lock.release()
            break
        data_name, data_x = expr_dict.popitem()
        data_data = data_x["data"]
        number_cells_rem = len(expr_dict)
        if number_cells_rem % 100 == 0:
            status_que.put({"print_out": str(number_cells_rem)})
            percent_complete = (int(data_x["numb_gene_set"]) - len(expr_dict)) / int(data_x["numb_gene_set"])*100
            percent_complete = round(percent_complete)
            print(percent_complete)
            status_que.put({"percent": percent_complete})
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


def producer(expr_dict, expr_obj, gmt_obj, lock, status_que):
    all_genes_list = matmulhypergeom.get_all_genes(expr_obj, gmt_obj)
    number_all_genes = len(all_genes_list)
    expmat, col_name, sample_expr_gene_num = matmulhypergeom.get_expmat(expr_obj, all_genes_list)
    gmtmat, row_name, pathway_gene_num = matmulhypergeom.get_gmtmat(gmt_obj, all_genes_list)
    intersect_matrix = matmulhypergeom.intersect_mat(expmat, gmtmat)
    # test_fewer = False
    # test_consumer = False
    for ip, rname in tqdm(enumerate(row_name)):
        _gene_set = pathway_gene_num[ip]
        expr_in_geneset_each_row = intersect_matrix[ip,]
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

    # print("update finished")
    status_que.put({"print_out": "update finished"})
    # consummer(expr_dict, gsea_samples, lock)
    # if test_consumer:
    #     consummer_row(expr_dict, gsea_samples, lock)


def pool_init(q, q_out):
    global queue  # make queue global in workers
    global queue_out
    queue = q
    queue_out = q_out


def show_w(para_que, status_que):
    app = QApplication(sys.argv)
    winf = main_win_setup.MyMainWindow(para_que=para_que, status_que=status_que)

    winf.show()
    sys.exit(app.exec_())


def receive_p(para_que_, status_que_, para_dic_):
    while True:
        try:
            cont = para_que_.get()
        except:
            raise ("queue error!")
        if len(cont) == 1 and para_que_.empty():
            if "finish" in cont:
                status_que_.put({"print_out": " get parameter finished\n start to run...\n"})
                break
        if len(cont) > 1:
            for k, v in cont.items():
                if k == "data_expr_dir":
                    para_dic_["batch_expr_dir"] = v
                if k == "gmt_dir":
                    para_dic_["batch_gmt_dir"] = v
                if k == "output_dir":
                    para_dic_["out_dir"] = v
                if k == "core_num":
                    para_dic_["core_num"] = v


if __name__ == '__main__':
    freeze_support()
    mgr = Manager()
    para_que = mgr.Queue()
    status_que = mgr.Queue()
    para_dic = mgr.dict()
    shw = Process(target=show_w, args=(para_que, status_que,))
    shw.start()
    recv = Process(target=receive_p, args=(para_que, status_que, para_dic,))
    recv.start()
    recv.join()

    batch_expr_dir = para_dic["batch_expr_dir"]
    batch_gmt_dir = para_dic["batch_gmt_dir"]
    out_dir = para_dic["out_dir"]
    n_core = para_dic["core_num"]

    for expr_each in os.listdir(batch_expr_dir):
        if os.path.isdir(os.path.join(batch_expr_dir, expr_each)):
            status_que.put({"print_out": "ready to processing data:"})
            expr_dir = os.path.join(batch_expr_dir, expr_each)
            status_que.put({"print_out": expr_dir+str("\n")})

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

                    lock = mgr.Lock()
                    expr_dict = mgr.dict()

                    gsea_samples = mgr.dict()

                    n_consumer = n_core
                    pool = Pool(n_consumer)

                    # print('ready to run multiprocessing')
                    # print('--' * 20)
                    # print('uploading data...')
                    status_que.put({"print_out": "ready to run multiprocessing"})
                    status_que.put({"print_out": '--' * 20})
                    status_que.put({"print_out": 'uploading data...'})

                    ti = time.time()
                    producer(expr_dict, expr_obj=expr_data, gmt_obj=gmt_data, lock=lock, status_que=status_que)
                    # prod = Process(target=producer, args=(expr_dict, expr_data, gmt_data,))
                    # prod.start()
                    # prod.join()
                    bg = len(expr_data.gene_list)
                    # print('run calculator...')
                    status_que.put({"print_out": "run calculator..."})

                    for i in range(n_consumer):
                        pool.apply_async(consummer_row, args=(expr_dict,
                                                              gsea_samples,
                                                              lock,
                                                              status_que))
                    status_que.put({"print_out": "starting:\n"})
                    # print("starting\n")
                    pool.close()
                    pool.join()
                    t1= time.time()
                    print("time:")
                    samples_gsea = matmulhypergeom.format_transfer_row(gsea_samples)
                    print("time:", time.time()-t1)
                    save_file(samples_gsea, formater, func_counts, sample_counts)
                    print("time:", time.time() - t1)
                    copyfile(os.path.join(expr_dir, 'barcodes.tsv'), os.path.join(each_out_dir, 'barcodes.tsv'))

                    pool.terminate()
                    # prod.terminate()
                    te = time.time()
                    # print("speed: {}cells/(core*sec)".format(len(expr_data.sample_geneset) / (te * n_consumer)))
                    status_que.put({"print_out": "speed: {}cells/(core*sec)".format(len(expr_data.sample_geneset) / ((te-ti) * n_consumer))})
                    del lock
                    del expr_dict
                    del samples_gsea
                    del expr_data
                    del gmt_data
                    del bg
                    del gsea_samples

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
                    status_que.put({"print_out": "time used:{}".format(te-ti)})
    status_que.put({"print_out": "all_task_finished!"})
    status_que.put({"end": "finished"})
    shw.join()
