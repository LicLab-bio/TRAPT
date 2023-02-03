# Acquire potential regulatory elements near that gene.
```shell
python3 geneMapper.py -i annotation/dhs_hg38_rose.bed -g HG38 -o output
```

# Calculating a regulatory potential matrix.
```shell
python3 calcRP.py -p ../../ -t chip -i output/dhs_hg38_rose_DHS_TO_GENE.txt -o output
python3 calcRP.py -p ../../ -t atac -i output/dhs_hg38_rose_DHS_TO_GENE.txt -o output
```
