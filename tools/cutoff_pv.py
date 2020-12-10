import math
import os
import sys

root_dir = "./"
data_file = "matrix.mtx"
tem_data_file = data_file + ".tem"


for each_dir in os.listdir(root_dir):
    sub_dir = os.path.join(root_dir, each_dir)
    print("current directory:", os.getcwd())
    os.chdir(sub_dir)
    with open(data_file, 'r') as sfp:
        with open(tem_data_file, 'w') as tfp:
            line_1 = sfp.readline()
            line_2 = sfp.readline()
            line_3_meta = sfp.readline()
            
            tfp.write(line_1)
            tfp.write(line_2)
            tfp.write(line_3_meta)
            
            for line in sfp.readlines():
                gene_idx, cell_idx, value = line.strip('\n').split(' ')
                value = float(value)
                value = round(value, 2)
                if value < 2:
                    value = 0.01
                tfp.write(' '.join([gene_idx, cell_idx, str(value)]) + '\n')
    os.remove(data_file)
    os.rename(tem_data_file, data_file)



