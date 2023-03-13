output=$1

if [ ! -f "$output/sample_names_chip.csv" ]
then
    echo "TSFS screening H3K27ac ChIP-seq samples."
    python3 feature_select.py $output chip
fi

if [ ! -f "$output/sample_names_atac.csv" ]
then
    echo "TSFS screening ATAC-seq samples."
    python3 feature_select.py $output atac
fi

if [ ! -f "$output/chip_dhs.npz" ]
then
    echo "Get H3K27ac-DHS matrix."
    python3 get_sample_dhs.py $output chip
fi

if [ ! -f "$output/atac_dhs.npz" ]
then
    echo "Get ATAC-DHS matrix."
    python3 get_sample_dhs.py $output atac
fi

if [ ! -f "$output/auc_tr_chip_ad.h5ad" ]
then
    echo "Get TR-H3K27ac AUC."
    python3 calc_auc_tr_chip.py $output
fi

if [ ! -f "$output/auc_tr_atac_ad.h5ad" ]
then
    echo "Get TR-ATAC AUC."
    python3 calc_auc_tr_atac.py $output
fi

if [ ! -f "$output/auc_tr_cre_ad.h5ad" ]
then
    echo "Get TR-CRE AUC."
    python3 calc_auc_tr_cre.py $output
fi

if [ ! -f "$output/auc_rp_chip_ad.h5ad" ]
then
    echo "Get RP-H3K27ac AUC."
    python3 calc_auc_rp_chip.py $output
fi

if [ ! -f "$output/auc_rp_atac_ad.h5ad" ]
then
    echo "Get RP-ATAC AUC."
    python3 calc_auc_rp_atac.py $output
fi

if [ ! -f "$output/STM.npz" ]
then
    echo "VGAE prediction TR interactive network."
    python3 calc_stm.py $output
fi

if [ ! -f "$output/chip.bw" ]
then
    echo "H3K27ac DHS bw file."
    python3 bed2bw.py $output chip
fi

if [ ! -f "$output/atac.bw" ]
then
    echo "ATAC DHS bw file."
    python3 bed2bw.py $output atac
fi

if [ ! -f "$output/input_gene.gff3.gz.tbi" ]
then
    echo "Input genes gff3 file."
    python3 gene2gff3.py $output
fi

if [ ! -f "$output/activity_summary.csv" ]
then
    echo "Estimated TR activity."
    python3 calc_activity.py $output
fi
