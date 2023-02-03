# # 获取共识性活性增强子(CAE)
import difflib
import os
import re
import sys

import pandas as pd

from args import args_

args = args_(sys.argv[1])

print(args.output)

def process_sample(sample_dict, sample):
    sample_ = re.sub(r'\._cas\.1$', '+', sample)
    sample_ = re.sub(r'\._cas$', '-', sample_)
    sample_ = re.sub(r'_cas$', '', sample_)
    sample_ = re.sub(r'^X', '', sample_)
    sample_ = re.sub(r'\.', '-', sample_)
    val = difflib.get_close_matches(sample_.upper(), list(sample_dict.keys()))
    if val:
        key = val[0]
        return sample,sample_dict.get(key)
    else:
        return sample,None

sample_bed = os.listdir('input/h3k27ac_bed_hg38')
sample_bed = pd.DataFrame(sample_bed,columns=['file_name'])
sample_bed['file_name_upper'] = sample_bed['file_name'].str.upper()

sample_dict = dict(sample_bed[['file_name_upper','file_name']].values.tolist())

samples = pd.read_csv('%s/sample_names_chip.csv' % args.output,header=None)[0]
samples = pd.DataFrame(map(lambda sample: process_sample(sample_dict, sample),samples.values),columns=['sample','file_name'])
samples = samples['file_name'].dropna().drop_duplicates()
data = []
for sample in samples:
    bed = pd.read_csv(f'input/h3k27ac_bed_hg38/{sample}',sep='\t',header=None)
    bed.index = [sample] * len(bed)
    data.append(bed)

data = pd.concat(data)
data = data.sort_values([0,1,2])
data.to_csv(f'{args.output}/samples.bed',header=False,index=False,sep='\t')
os.system(f'bedtools merge -c 4 -o count -i {args.output}/samples.bed > {args.output}/samples_merge.bed;' \
    f'rm {args.output}/samples.bed')
samples_merge = pd.read_csv(f'{args.output}/samples_merge.bed',sep='\t',header=None)
CAE = samples_merge[samples_merge[3] > len(samples) * 0.5]
CAE.index = list(map(lambda i:f'CAE_{i}',range(len(CAE))))
CAE.to_csv(f'{args.output}/CAE.bed',header=False,sep='\t')
