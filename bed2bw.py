# http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
# http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
import os
import re
from sys import argv

import pandas as pd
from scipy import sparse

from args import args_

# python3 bed2bb.py chip AHR@DataSet_01_283_up500


# /home/tostring/rgtdata/hg38/genome_hg38.fa
# /data/zgr/data/bulk数据处理/ref/hg19.fasta
output = argv[1]
args = args_(output)
sample = argv[2]
# /home/tostring/.local/bin/CrossMap.py
# crossmap=argv[2]

# significance_index = pd.read_csv(f'{output}/significance_index_{sample}.csv',header=None)[0].values

sample_dhs = sparse.load_npz(f'{output}/{sample}_dhs_pred.npz')
sample_dhs = pd.DataFrame(sample_dhs.T.toarray())

dhs = pd.read_csv('input/dhs_hg38.bed',sep='\t',header=None) 
dhs = dhs.iloc[:,:3].reset_index().drop(['index'],axis=1)

# 坐标转换
# pip install CrossMap
# wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -O input/hg19ToHg38.over.chain.gz

# os.system(f"""
# python3 {crossmap} bed input/hg19ToHg38.over.chain.gz input/dhs.bed input/dhs_hg38.bed
# """)

bedgraph = pd.concat([dhs,sample_dhs],axis=1)
bedgraph.columns = [0,1,2,3]
bedgraph = bedgraph[bedgraph[2]>bedgraph[1]]
bedgraph = bedgraph.sort_values(by=[0,1,2],ignore_index=True).astype(str)
bedgraph = bedgraph[bedgraph[0].apply(lambda x:bool(re.match('^chr[\dX-Y]{1,2}$',x)))]

print(bedgraph)
bedgraph.to_csv(f'{output}/{sample}.bed',index=False,header=False,sep='\t')

os.system(f"""
# 区域去重叠
bedtools merge -i {output}/{sample}.bed > {output}/{sample}_merge.bed
bedtools intersect -a {output}/{sample}_merge.bed -b {output}/{sample}.bed \\
-wa -wb|bedtools groupby -i - -g 1-3 -c 7 -o sum > {output}/{sample}.bedgraph
rm -f {output}/{{{sample},{sample}}}_merge.bed

awk '{{if ($4 > 1) $4=1}}1' {output}/{sample}.bedgraph > {output}/{sample}_modify.bedgraph
bedGraphToBigWig {output}/{sample}_modify.bedgraph input/hg38.chrom.sizes {output}/{sample}.bw
# rm -f {output}/{{{sample},{sample}}}_modify.bedgraph
""")
