import sys
import os
import numpy as np
import pandas as pd
from abc import ABCMeta, abstractmethod


sys.path.append('./')

def read_matrix_data(filename, sep='\t'):
    pass

class MatrixObj():
    def __init__(self, filename, sep='\t'):
        self.fn = filename
        self.sep = sep

    @abstractmethod
    def read_data(self):
        pass

    @abstractmethod
    def write_10x(self):
        pass

if __name__ == '__main__':
    pass




'''
。                         鸡你太美
                          鸡你实在太美
                           鸡你是太美
                              鸡你太美
                         实在是太美鸡你
               鸡你 实在是太美鸡你  美
            鸡你  实在是太美鸡美   太美
        鸡你     实在是太美鸡美      太美        
    鸡你       实在是太美鸡美          太美
  鸡你       鸡你实在是美太美       美蓝球球
鸡 鸡       鸡你实在是太美          篮球篮球球
鸡            鸡你太美裆鸡太啊         蓝篮球
                鸡你太美裆裆鸡美
                  鸡你美裆    裆鸡美
                   鸡太美         鸡太美
                     鸡美              鸡美
                     鸡美                鸡美
                      鸡美                鸡美
                      鸡太                鸡太
                    金 猴                金猴
                    皮 鞋                皮鞋金猴
                     金光                金光 大道
                    大道

'''
