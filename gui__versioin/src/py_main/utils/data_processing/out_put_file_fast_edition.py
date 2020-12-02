import math
import shutil
import os
import sys


class Formater():
    def __init__(self, meta_data,  out_dir, funclist, thredhold=0.5):
        #original matrix meta data
        self.meta_data = meta_data
        self.thredhold = thredhold

        #gsea data
        self.funclist = funclist

        #out put
        self.fp_matrix = os.path.join(out_dir, "matrix.mtx")
        self.fp_functions = os.path.join(out_dir, "genes.tsv")

        #save record counts
        self.record_counts = 0

    def to_orig_marmat_matrix_meta(self, func_counts = 0, sample_counts = 0):
        with open(self.fp_matrix, 'r') as fp:
            contents = fp.read()

        cache = ""
        with open(self.fp_matrix, 'w') as wf:
            cache += '\n'.join([self.meta_data[0], self.meta_data[1]]) + '\n'
            wf.write(cache)
            spn = " "
            metas = spn.join([str(func_counts), str(sample_counts), str(self.record_counts)])
            wf.write(metas)
            wf.write('\n')
            wf.write(contents)

    def to_orig_format_matrix_new(self, sample_gsea):
        self.matrix_cache = ""
        func_dic = {k: v for v, k in enumerate(self.funclist) }

        with open(self.fp_matrix, 'a') as wf:
            for s_id, each_p in sample_gsea.items():
                for gene_set_name, pv in each_p:
                    pv = float(pv)
                    if (pv >= self.thredhold):
                        continue
                    pv = -math.log(pv, 10)
                    term_idx = func_dic[gene_set_name]
                    # term_idx = self.funclist.index(gene_set_name)
                    cache += " ".join([str(term_idx + 1), str(s_id + 1), str(pv)]) + "\n"
                    self.record_counts += 1

            wf.write(cache)
    def processing_to_orig_format_matrix_new(self, sample_gsea):
        pass

    def to_orig_format_func(self):
        cache = ""
        with open(self.fp_functions, 'w') as wf:
            for f_item in self.funclist:
                cache += "\t".join(["function_items", f_item]) + '\n'
            wf.write(cache)

    def to_matrix_format(self):
        pass







