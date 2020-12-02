import os
from collections import OrderedDict

class GetData():
    def __init__(self, data_dir):
        self.dir = data_dir
        self.gene_list = None
        self.matrix_ex = None
        self.matrix_ex_meta_data = None
        self.matrix_ex_meta_data_head = None
        self.sample_geneset = None

    # def get_bacode(self,):
    #     pass
    def get_genename(self, _genes):
        gene_list = []

        spline = _genes[0].split("\t")
        if len(spline) == 1:
            for gene in _genes:
                gene_name = gene
                gene_list.append(gene_name)
            self.gene_list = gene_list
        elif len(spline) == 2:
            for gene in _genes:
                _i, gene_name = gene.split("\t")
                gene_list.append(gene_name)
            self.gene_list = gene_list
        else:
            raise ("gene file format can't recognize!")

    def get_matrix(self, _matrix):
        meta_data_head = _matrix[0:2]
        del _matrix[0:2]
        meta_data = _matrix[0]
        n_gene, n_sample, n_items = meta_data.split(" ")
        meta_data = (n_gene, n_sample, n_items)
        self.matrix_ex_meta_data = meta_data
        self.matrix_ex_meta_data_head = meta_data_head
        del _matrix[0]
        matrix_ex = []
        for line in _matrix:
            gene_idex, sample_id, value = line.split(" ")

            gene_idex = int(gene_idex)
            sample_id = int(sample_id)
            value = int(round(float(value)))

            gene_idex = gene_idex - 1
            sample_id = sample_id - 1
            matrix_ex.append((gene_idex, sample_id, value))
        self.matrix_ex = matrix_ex

    def read_data(self, file_name, _lower=False):
        fn = file_name
        with open(fn, 'r') as fp:
            file_lines = []
            for line in fp.readlines():
                if _lower:
                    line = line.lower()

                line = line.strip('\n')
                file_lines.append(line)
        return file_lines

    def get_sample2geneset(self):
        sample_geneset = OrderedDict()
        for i in range(int(self.matrix_ex_meta_data[1])):
            sample_geneset[i] = set()
        for tp in self.matrix_ex:
            gene_idex, sample_id, _ = tp
            sample_geneset[sample_id].add(self.gene_list[int(gene_idex)])
        self.sample_geneset = sample_geneset

    def run_1(self):
        gene_file = "genes.tsv"
        mtx_file = "matrix.mtx"
        gene_file = os.path.join(self.dir, gene_file)
        mtx_file = os.path.join(self.dir, mtx_file)

        fl = self.read_data(gene_file)
        self.get_genename(fl)

        fl = self.read_data(mtx_file)
        self.get_matrix(fl)

        self.get_sample2geneset()


class Gmt_stat():
    def __init__(self, file_name, _lower=True):
        self.file_name = file_name
        self.term_set = None
        self.term_quantity = None
        self.term_names = None
        self.term_all_gene = None
        self.background_len = None
        self._lower = _lower
        self.term_number = {}
        self.gene_term_dic = {}

    def get_gmtobj(self):
        term_set = OrderedDict()
        term_names = []
        term_all_gene = set()
        with open(self.file_name, 'r') as fgmt:
            for line in fgmt.readlines():
                line = line.strip('\n')
                if self._lower:
                    line = line.lower()
                term_name, _, *gene_set = line.split('\t')
                term_name = term_name.lower()
                term_set[term_name] = set(gene_set)
                term_all_gene.update(set(gene_set))
                term_names.append(term_name)
        self.term_set = term_set
        self.term_names = term_names
        self.term_all_gene = term_all_gene

    def get_gmt_info(self):
        self.term_quantity = len(self.term_set)
        self.background_len = len(self.term_all_gene)
        #self.get_dict()

    def get_dict(self):
        for term, _geneset in self.term_set.items():
            self.term_number[term] = len(_geneset)
            for gene in _geneset:
                if gene in self.gene_term_dic:
                    self.gene_term_dic[gene].add(term)
                else:
                    self.gene_term_dic[gene] = set(gene)
