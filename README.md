### TRAPT

[TRAPT](http://www.licpathway.net/TRAPT) is a novel deep learning framework for transcription regulators prediction via integraing large-scale epigenomic data.

### Usage

First, download library: 

[TRAPT library](http://www.licpathway.net/TRAPT/download)


Second, install TRAPT:

```sh
conda create --name TRAPT python=3.7
conda activate TRAPT
pip install TRAPT
```

Run TRAPT using a [case](http://www.licpathway.net/TRAPT/static/download/ESR1@DataSet_01_111_down500.txt):

Use shell commands:
```bash
# help
trapt --help
# run
trapt --library library \
      --input ESR1@DataSet_01_111_down500.txt \
      --output output/test/ESR1@DataSet_01_111_down500
```

Using the python interface:
```python
import os
from TRAPT.Tools import Args, RP_Matrix
from TRAPT.Run import runTRAPT

# library path
library = 'library'
# input file path
input = 'ESR1@DataSet_01_111_down500.txt'
# output file path
output = 'output/test/ESR1@DataSet_01_111_down500'

rp_matrix = RP_Matrix(library)
args = Args(input, output)
runTRAPT([rp_matrix, args])
```

### De novo library

```bash
# Constructing TR-RP matrix
python3 CalcTRRPMatrix.py library
# Constructing H3K27ac-RP matrix
python3 CalcSampleRPMatrix.py H3K27ac library
# Constructing ATAC-RP matrix
python3 CalcSampleRPMatrix.py ATAC library
# Reconstruct TR-H3K27ac adjacency matrix
python3 DLVGAE.py H3K27ac library
# Reconstruct TR-ATAC adjacency matrix
python3 DLVGAE.py ATAC library
# Prediction (TR-H3K27ac)-RP matrix
python3 CalcTRSampleRPMatrix.py H3K27ac library
# Prediction (TR-ATAC)-RP matrix
python3 CalcTRSampleRPMatrix.py ATAC library
```