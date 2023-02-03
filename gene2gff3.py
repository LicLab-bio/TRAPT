import csv
import os
import re
from sys import argv

import pandas as pd

from args import args_

output = argv[1]
args = args_(output)
gene_set = set(pd.read_csv(args.input,header=None)[0])
gencode_gene = pd.read_csv('input/gencode.v40.annotation.sorted.gff3.gz',header=None,sep="|",on_bad_lines='error')
genes = gencode_gene[0].apply(lambda x:''.join(re.findall(r'gene_name=(.*?);',x)))
genes = genes.apply(lambda x: x in gene_set)
input_gene = gencode_gene[genes]
input_gene.to_csv(f'{output}/input_gene.gff3',sep='|',header=False,index=False,quoting=csv.QUOTE_NONE,quotechar="|")
os.system(f"""
bgzip -f {output}/input_gene.gff3
tabix -f {output}/input_gene.gff3.gz
""")
