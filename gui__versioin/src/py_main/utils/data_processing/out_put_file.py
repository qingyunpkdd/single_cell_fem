import math
import shutil
import os

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

        with open(self.fp_matrix, 'w') as fp:
            fp.write(self.meta_data[0])
            fp.write('\n')
            fp.write(self.meta_data[1])
            fp.write('\n')
            #original file meta information
            # spn = " "
            # mt = spn.join(list(self.data.matrix_ex_meta_data))

            spn = " "
            metas = spn.join([str(func_counts), str(sample_counts), str(self.record_counts)])
            fp.write(metas)
            fp.write('\n')
            fp.write(contents)

    def to_orig_format_matrix(self,sample_gsea):

        with open(self.fp_matrix, 'a') as fp:

            for s_id, dfs in sample_gsea.items():
                for index, row in dfs.iterrows():
                    pv = row["Adjusted P-value"]
                    pv = float(pv)
                    #                    print(pv)
                    if (pv >= self.thredhold):
                        continue
                    pv = -math.log(pv, 10)
                    #                    print("-log_pvalue:",pv)
                    term = row["Term"]
                    term_idx = self.funclist.index(term)
                    fp.write(str(term_idx + 1))
                    fp.write(" ")
                    fp.write(str(s_id + 1))
                    fp.write(" ")
                    fp.write(str(pv))
                    fp.write('\n')
                    self.record_counts += 1
    def to_orig_format_matrix_new(self, sample_gsea):
        func_dic = {k: v for v, k in enumerate(self.funclist) }
        with open(self.fp_matrix, 'a') as fp:
            for s_id, each_p in sample_gsea.items():
                for gene_set_name, pv in each_p:
                    pv = float(pv)
                    if (pv >= self.thredhold):
                        continue
                    pv = -math.log(pv, 10)
                    term_idx = func_dic[gene_set_name]
                    fp.write(str(term_idx + 1))
                    fp.write(" ")
                    fp.write(str(s_id + 1))
                    fp.write(" ")
                    fp.write(str(pv))
                    fp.write('\n')
                    self.record_counts += 1

    def to_orig_format_func(self):
        with open(self.fp_functions, 'w') as fp:
            for f_item in self.funclist:
                fp.write("function_items")
                fp.write("\t")
                fp.write(f_item)
                fp.write('\n')

    def to_matrix_format(self):
        pass