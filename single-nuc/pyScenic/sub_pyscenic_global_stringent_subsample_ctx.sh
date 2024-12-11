#PBS -lselect=1:ncpus=100:mem=980gb
#PBS -lwalltime=48:00:00

export PATH=$Home/rds/general/user/username/home/.conda/envs/personal/bin/python:$PATH
cd path_in

module load anaconda3/personal
unset PYTHONPATH
source activate py_sc2023

pyscenic ctx /path_in/adj_merged.tsv \ #obtained from pyscenic_create_mean_adjacency.py
    encode_20190621__ChIP_seq_transcription_factor.hg38__refseq-r80__10kb_up_and_down_tss.max.genes_vs_tracks.rankings.feather \
    --annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname global_all_stringent.loom \ # from pyscenic_data_prep_stringent
    --output  /path_out/reg_global_stringent_subsample.csv \ 
    --mask_dropouts \
    --min_genes 5 \
    --no_pruning \
    --chunk_size 50 \
    --all_modules \
    --num_workers 20