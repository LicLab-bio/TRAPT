import json
import os
import subprocess

import pandas as pd


class args_:
    def __init__(self,path) -> None:
        self.parmas_json = os.path.join(path,'params.json')
        self.data = self.get_params()
        self.input = self.data['input']
        self.output = self.data['output']
        self.dim = 128
        self.top_n = 1
        self.min_inter_limit = 20
        self.vgae_iterations = 10
        self.range_threshold = 0.05
        if 'vgae_iterations' in self.data:
            self.vgae_iterations = self.data['vgae_iterations']
        if 'range_threshold' in self.data:
            self.range_threshold = self.data['range_threshold']
    
    @staticmethod
    def get_wc():
        wc = subprocess.run("ps -ef|grep -E 'calc_stm.py|_dhs.py|feature_select'|wc -l", shell=True,capture_output=True)
        # print(subprocess.run("ps -ef|grep -E 'calc_stm.py|_dhs.py|feature_select'",shell=True,capture_output=True).stdout.decode('utf8'))
        wc = int(wc.stdout.decode('utf8').replace('\n','')) - 3
        return wc

    def get_params(self):
        with open(self.parmas_json,'r') as file:
            data = json.loads(file.read())
        return data

    def get_dhs_index(self):
        input_genes = pd.read_csv(self.input,header=None)
        input_genes.columns = ['gene']
        dhs_gene = pd.read_csv('input/dhs_gene_hg38.csv',sep='\t')
        # Only DHS near that set of us input genes are consider
        input_genes_dhs = input_genes.merge(dhs_gene,on='gene')
        return input_genes_dhs['index'].drop_duplicates().values
