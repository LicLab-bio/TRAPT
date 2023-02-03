import functools
import glob
import json
import os
import re
from uuid import uuid1

import pandas as pd

PATH = "/mnt/data2/python_web/TCAPT"
PATH_OUT = f"{PATH}/jbrowse2/output"

def mException(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            print(e)
            return None

    return wrapper


def run(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        os.chdir(PATH)
        try:
            return func(*args, **kwargs)
        except Exception as e:
            print(e)
            return None
        finally:
            os.chdir(PATH_OUT)

    return wrapper

@mException
def _getTRsActivityTable(*args, **kwargs):
    uuid = kwargs['uuid']
    data = pd.read_csv(f'{PATH_OUT}/{uuid}/activity_summary.csv')
    data = data.round(4).reset_index()
    data.columns = ['Rank','TR name','TR activity','Interaction score','CRE score','ATAC score','H3K27ac score','RP ATAC score','RP H3K27ac score']
    data['Rank'] = data['Rank'] + 1
    r = data.dropna().values.tolist()
    return r

@mException
def _isFinal(*args, **kwargs):
    uuid = kwargs['uuid']
    is_final = os.path.exists(f'{PATH_OUT}/{uuid}/activity_summary.csv')
    return is_final

def _drawLine(*args, **kwargs):
    uuid = kwargs['uuid']
    tr = kwargs['tr']
    views_tracks = ['input_gene','predicted_chip_signal','predicted_atac_signal']
    with open(f'{PATH_OUT}/{uuid}/tr_top_samples.json','r') as file:
        tr_top_samples = json.loads(file.read())
        views_tracks.extend(tr_top_samples.get(tr)[:5])
    tr_all = glob.glob(f'{PATH}/jbrowse2/data/TRs/{tr}@*.sorted.bed.gz')
    tr_tracks = list(map(lambda x:re.findall('TRs/(.*?)\.sorted',x)[0],tr_all))
    # views_tracks.extend(tr_tracks[:5])

    tracks = []
    tracks.extend([
        {
            "type": "FeatureTrack",
            "trackId": "input_gene",
            "name": "input_gene",
            "category": ["Annotation"],
            "assemblyNames": ["genome_hg38"],
            "adapter": {
                "type": "Gff3TabixAdapter",
                "gffGzLocation": {
                "uri": "input_gene.gff3.gz",
                "locationType": "UriLocation"
                },
                "index": {
                "location": {
                    "uri": "input_gene.gff3.gz.tbi",
                    "locationType": "UriLocation"
                }
                }
            },
        },
        {
            "type": "FeatureTrack",
            "trackId": "gencode_gene",
            "name": "gencode_gene",
            "category": ["Annotation"],
            "assemblyNames": ["genome_hg38"],
            "adapter": {
                "type": "Gff3TabixAdapter",
                "gffGzLocation": {
                "uri": f"../../data/gencode.v40.annotation.sorted.gff3.gz",
                "locationType": "UriLocation"
                },
                "index": {
                "location": {
                    "uri": f"../../data/gencode.v40.annotation.sorted.gff3.gz.tbi",
                    "locationType": "UriLocation"
                }
                }
            },
        },
        {
            "type": "FeatureTrack",
            "trackId": "DHS_hg38",
            "name": "DHS_hg38",
            "category": ["Annotation"],
            "assemblyNames": ["genome_hg38"],
            "adapter": {
                "type": "BedTabixAdapter",
                "bedGzLocation": {
                "uri": f"../../data/dhs_hg38.sorted.bed.gz",
                "locationType": "UriLocation"
                },
                "index": {
                "location": {
                    "uri": f"../../data/dhs_hg38.sorted.bed.gz.tbi",
                    "locationType": "UriLocation"
                },
                "indexType": "TBI"
                }
            },
        },
    ])
    tracks.extend([
        {
        "type": "QuantitativeTrack",
        "trackId": "predicted_chip_signal",
        "name": "predicted_chip_signal",
        "category": ["predicted_signal_BigWig"],
        "adapter": {
            "type": "BigWigAdapter",
            "bigWigLocation": {
            "uri": f"chip.bw",
            "locationType": "UriLocation"
            }
        },
        "assemblyNames": [
            "genome_hg38"
        ]
        },
        {
        "type": "QuantitativeTrack",
        "trackId": "predicted_atac_signal",
        "name": "predicted_atac_signal",
        "category": ["predicted_signal_BigWig"],
        "adapter": {
            "type": "BigWigAdapter",
            "bigWigLocation": {
            "uri": f"atac.bw",
            "locationType": "UriLocation"
            }
        },
        "assemblyNames": [
            "genome_hg38"
        ]
        },
    ])
    for tr in tr_tracks:
        tracks.append({
        "type": "QuantitativeTrack",
        "trackId": tr,
        "name": tr,
        "category": ["TRs BigWig"],
        "adapter": {
            "type": "BigWigAdapter",
            "bigWigLocation": {
            "uri": f"../../data/TRs/{tr}.bw",
            "locationType": "UriLocation"
            }
        },
        "assemblyNames": [
            "genome_hg38"
        ]
        })

    config_json = {
    "assemblies": [
    {
        "name": "genome_hg38",
        "sequence": {
        "type": "ReferenceSequenceTrack",
        "trackId": "genome_hg38",
        "adapter": {
            "type": "IndexedFastaAdapter",
            "fastaLocation": {
            "uri": f"../../data/genome_hg38.fa",
            "locationType": "UriLocation"
            },
            "faiLocation": {
            "uri": f"../../data/genome_hg38.fa.fai",
            "locationType": "UriLocation"
            }
        }
        }
    }
    ],
    "tracks": tracks
    }
    params = {
        "views": [
            {
                "type": "LinearGenomeView",
                "loc": "chr1:1-248956422",
                "assembly": "genome_hg38",
                "tracks": views_tracks,
            }
        ]
    }
    with open(f'{PATH_OUT}/{uuid}/config.json','w') as file:
        file.write(json.dumps(config_json))
    return params

@run
def _run(*args, **kwargs):
    uuid = uuid1()
    vgae_iterations = kwargs['vgae_iterations']
    range_threshold = kwargs['range_threshold']
    gene_set = re.sub(r'[\s,;]+','\n',kwargs['gene_set'])
    os.system(f'mkdir -p {PATH_OUT}/{uuid}')
    with open(f'{PATH_OUT}/{uuid}/gene_set.txt','w') as file:
        file.write(gene_set)
    with open(f'{PATH_OUT}/{uuid}/params.json','w') as file:
        params = {
            'vgae_iterations': int(vgae_iterations),
            'range_threshold': float(range_threshold),
            'input': f'{PATH_OUT}/{uuid}/gene_set.txt',
            'output': f'{PATH_OUT}/{uuid}'
        }
        file.write(json.dumps(params))
    os.system(f'sh run.sh {PATH_OUT}/{uuid} > {PATH_OUT}/{uuid}/run.log 2>&1 &')
    print(uuid)
    return uuid

@run
def _tryAgain(*args, **kwargs):
    uuid = kwargs['uuid']
    os.system(f'sh run.sh {PATH_OUT}/{uuid} > {PATH_OUT}/{uuid}/run.log 2>&1 &')
    return True

@mException
def _getLog(*args, **kwargs):
    uuid = kwargs['uuid']
    with open(f'{PATH_OUT}/{uuid}/run.log','r') as file:
        log = file.read()
        log = re.sub('\n','<br>',log)
        return log
