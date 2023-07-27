import functools
import os

import numpy as np
import pandas as pd

from TRAPT.Tools import RP_Matrix, Args, Type
from TRAPT.CalcTRAUC import CalcTRAUC
from TRAPT.DLFS import FeatureSelection

def get_params(func):
    @functools.wraps(func)
    def wrapper(args):
        return func(*args)

    return wrapper


@get_params
def runTRAPT(rp_matrix: RP_Matrix, args: Args):
    obs = rp_matrix.TR.obs

    if os.path.exists(f"{args.output}/H3K27ac_RP.csv"):
        H3K27ac_RP = pd.read_csv(f"{args.output}/H3K27ac_RP.csv", header=None)[0]
    else:
        FS_H3K27ac = FeatureSelection(args, rp_matrix.H3K27ac, Type.H3K27ac)
        H3K27ac_RP = FS_H3K27ac.run()
        H3K27ac_RP.to_csv(f"{args.output}/H3K27ac_RP.csv", index=False, header=False)

    if os.path.exists(f"{args.output}/ATAC_RP.csv"):
        ATAC_RP = pd.read_csv(f"{args.output}/ATAC_RP.csv", header=None)[0]
    else:
        FS_ATAC = FeatureSelection(args, rp_matrix.ATAC, Type.ATAC)
        ATAC_RP = FS_ATAC.run()
        ATAC_RP.to_csv(f"{args.output}/ATAC_RP.csv", index=False, header=False)

    if os.path.exists(f"{args.output}/RP_TR_H3K27ac_auc.csv"):
        RP_TR_H3K27ac_auc = pd.read_csv(
            f"{args.output}/RP_TR_H3K27ac_auc.csv", index_col=0, header=None
        )
    else:
        H3K27ac_RP = H3K27ac_RP.values.flatten()
        CTR_TR = CalcTRAUC(args, rp_matrix.TR_H3K27ac, H3K27ac_RP)
        RP_TR_H3K27ac_auc = CTR_TR.run()
        RP_TR_H3K27ac_auc.to_csv(f"{args.output}/RP_TR_H3K27ac_auc.csv", header=False)

    if os.path.exists(f"{args.output}/RP_TR_ATAC_auc.csv"):
        RP_TR_ATAC_auc = pd.read_csv(
            f"{args.output}/RP_TR_ATAC_auc.csv", index_col=0, header=None
        )
    else:
        ATAC_RP = ATAC_RP.values.flatten()
        CTR_TR = CalcTRAUC(args, rp_matrix.TR_ATAC, ATAC_RP)
        RP_TR_ATAC_auc = CTR_TR.run()
        RP_TR_ATAC_auc.to_csv(f"{args.output}/RP_TR_ATAC_auc.csv", header=False)

    data_auc = pd.concat([RP_TR_H3K27ac_auc, RP_TR_ATAC_auc], axis=1)
    data_auc /= np.linalg.norm(data_auc, axis=0, keepdims=True)
    TR_activity = pd.DataFrame(
        np.sum(data_auc.values, axis=1), index=data_auc.index, columns=[1]
    )
    TR_detail = pd.concat([TR_activity, data_auc], axis=1).reset_index()
    TR_detail.columns = ["TR", "TR activity", "RP_TR_H3K27ac_auc", "RP_TR_ATAC_auc"]
    obs.index.name = "TR"
    TR_detail = TR_detail.merge(obs.reset_index(), on="TR").sort_values(
        "TR activity", ascending=False
    )
    TR_detail.to_csv(os.path.join(args.output, "TR_detail.txt"), index=False, sep="\t")
    return TR_detail


def main():
    from optparse import OptionParser
    usage = "usage: %prog [options] -l [LIBRARY] -i [INPUT] -o [OUTPUT]"
    parser = OptionParser(usage = usage)
    parser.add_option("-l","--library", dest="library",nargs = 1, default="library", help = "Enter the library path, default is './library'")
    parser.add_option("-p","--threads", dest="threads",nargs = 1, default=16, help = "Number of threads to launch, default is 16")
    parser.add_option("-t","--trunk_size", dest="trunk_size",nargs = 1, default=32768, help = "Block size. If the memory is insufficient, set a smaller value. The default value is 32768")
    parser.add_option("-i","--input", dest="input",nargs = 1, default=None, help = "Enter a gene list")
    parser.add_option("-o","--output", dest="output",nargs = 1, default=None, help = "Enter an output folder")
    options,args = parser.parse_args()
    library = options.library
    threads = options.threads
    trunk_size = options.trunk_size
    input = options.input
    output = options.output
    rp_matrix = RP_Matrix(library)
    args = Args(input, output, threads, trunk_size)
    runTRAPT([rp_matrix, args])

if __name__ == "__main__":
    main()
