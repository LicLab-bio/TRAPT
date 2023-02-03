import re
from glob import glob
from multiprocessing import Pool
from optparse import OptionParser

import anndata
import numpy as np
import pandas as pd
from scipy import sparse
from tqdm import tqdm

usage = "usage: %prog [options] -p [PATH] -t [TYPE] -i [DHS_TO_GENE_FILE] -o [OUTPUT]"
parser = OptionParser(usage = usage)
parser.add_option("-p", dest="PATH",nargs = 1, default=None,help = "TCAPT tool path.")
parser.add_option("-t", dest="TYPE",nargs = 1, default=None,help = "(atac,chip)")
parser.add_option("-i", dest="DHS_TO_GENE_FILE",nargs = 1, default=None,help = "dhs_hg38_rose_DHS_TO_GENE.txt")
parser.add_option("-o", dest="OUTPUT",nargs = 1, default=None,help = "Enter an output folder.")
options,args = parser.parse_args()

ucsc_gene = pd.read_csv('annotation/hg38_refseq.ucsc',sep='\t')
ucsc_gene = ucsc_gene[['name','name2']]
ucsc_gene.columns = ['GENES','SYMBOLS']
ucsc_gene = ucsc_gene.drop_duplicates()
dhs_hg38 = pd.read_csv('annotation/dhs_hg38_rose.bed',sep='\t',header=None)
dhs_hg38 = dhs_hg38.reset_index()[['index',0]]
dhs_hg38.columns = ['index','DHS']
dhs_gene = pd.read_csv(options.DHS_TO_GENE_FILE,sep='\t')
dhs_gene = dhs_gene.iloc[:,[0,4,5]]
dhs_gene.columns = ['DHS','GENES','DISTANCE']
dhs_gene_merge = dhs_gene.merge(dhs_hg38,on='DHS')
dhs_gene_merge = dhs_gene_merge.merge(ucsc_gene,on='GENES')
dhs_gene_merge = dhs_gene_merge.drop_duplicates(['DHS','SYMBOLS'])
dhs_gene_g = dhs_gene_merge.groupby('SYMBOLS').agg({
    'DISTANCE':list,
    'index': list
})
dhs_gene_g['DISTANCE'] = dhs_gene_g['DISTANCE'].apply(lambda x:np.array(x))
genes = dhs_gene_g.index

def dhs2gene(args):
    sample,alpha = args
    try:
        sample_name = re.findall(r'dhs@(.*?).csv',sample)[0]
        vec = pd.read_csv(sample,sep='\t',header=None)[0].values
        rp_vec = []
        for gene in genes:
            index = dhs_gene_g.loc[gene,'index']
            s = vec[index]
            d = dhs_gene_g.loc[gene,'DISTANCE']
            w = (np.exp(-d/alpha) + 1)/(2 * np.exp(d/alpha))
            rp = np.mean(np.multiply(w,s))
            rp_vec.append(rp)
        rp_vec = np.log(np.array(rp_vec) + 1)
        rp_vec = rp_vec - np.median(rp_vec)
        rp_vec = sparse.csr_matrix(rp_vec)
        return sample_name,rp_vec
    except:
        print('Error %s !' % sample)
        return None

regu_matrix = []
samples = glob(f'{options.PATH}/input/{options.TYPE}_matrix/dhs@*_hg38.csv')

args = list(map(lambda sample:(sample,100000),samples))
regu_matrix = []
sample_list = []
with Pool(24) as pool:
    for row in tqdm(pool.imap(dhs2gene, args),total=len(args)):
        if row:
            sample_name,rp_vec = row
            sample_list.append(sample_name)
            regu_matrix.append(rp_vec)


regu_matrix = sparse.vstack(regu_matrix,dtype='float32')
regu_matrix_ad = anndata.AnnData(regu_matrix)
regu_matrix_ad.var_names = genes
regu_matrix_ad.obs_names = sample_list
regu_matrix_ad.write_h5ad(f'{options.OUTPUT}/regu_matrix_{options.TYPE}_ad.h5ad')
