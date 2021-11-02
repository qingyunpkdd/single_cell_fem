import os
# import sys
# import click
from tqdm import tqdm

# @click.command()
#
# @click.option('--file_dir', prompt='please input the input file', help='the file to be process as 10x format!')
# @click.option('--save_dir', prompt='please input the output dir', help='the output directory!')

class ToTenx(object):
    def __init__(self, file_dir, save_dir, f_sepr='\t',  matrix_head='%%MatrixMarket matrix coordinate real general\n%\n'):
        super().__init__()
        self.fdir = file_dir
        self.save_dir = save_dir
        self.matrix_head = matrix_head
        self.matrix_head_sta = None

        self.f_sepr = f_sepr

        self.barcodes = []        
        self.genelist = []
        self.matrix = []
        self.namedic = {'barcode':'barcodes.tsv', 'gene':'genes.tsv', 'matrix':'matrix.mtx'}
        self.matrx10x = []

        self.read()
        self.to10x()
        self.save()
        
    def read(self):
        with open(self.fdir,'r') as fp:
            print('read file...')
            #read the cell barcodes
            barcodes = fp.readline().strip('\n').split(self.f_sepr)
            self.barcodes = barcodes
            for line in tqdm(fp.readlines()):
                genename, *values = line.strip('\n').split(self.f_sepr)
                self.genelist.append(genename)
                self.matrix.append(values)
    
    def to10x(self):
        print('trim the zero!')
        for idi, gene_allcell in tqdm(enumerate(self.matrix)):
            for idj, value in enumerate(gene_allcell):
                if value == '0':
                    continue
                self.matrx10x.append([str(idi+1), str(idj+1), value])
        self.matrix_head_sta = [str(len(self.genelist)), str(len(self.barcodes)), str(len(self.matrx10x))]
        print('processed!')

    def save(self):
        #save barcode
        print('write to 10x')
        with open(os.path.join(self.save_dir, self.namedic['barcode']), 'w') as bf:
            content = '\n'.join(self.barcodes) + '\n'
            bf.write(content)
        with open(os.path.join(self.save_dir, self.namedic['gene']), 'w') as gf:
            content = None
            content = '\n'.join(self.genelist) + '\n'
            gf.write(content)
        with open(os.path.join(self.save_dir, self.namedic['matrix']), 'w') as mf:
            content = None
            mf.write(self.matrix_head)
            mf.write(' '.join(self.matrix_head_sta) + '\n')
            for line in tqdm(self.matrx10x):
                line = ' '.join(line)
                mf.write(line)
                mf.write('\n')

if __name__ == '__main__':
    inputfile = './totenxtest/test.txt'
    outputdir = './totenxtest/testout'

    totenx = ToTenx(inputfile, outputdir)
    #inputfile = inputfile
    #output = output


    #totenx = ToTenx('./totenx_test/test.txt', './totenx_test/testout/')



