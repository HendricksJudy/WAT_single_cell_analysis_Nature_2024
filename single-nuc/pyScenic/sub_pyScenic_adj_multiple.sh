#!/usr/bin/bash

#PBS -lselect=1:ncpus=200:mem=4000gb
#PBS -lwalltime=72:00:00

export PATH=$Home/rds/general/user/sername/home/.conda/envs/personal/bin/python:$PATH
cd path_in
module load anaconda3/personal
unset PYTHONPATH
source activate py_sc2023

arboreto_with_multiprocessing.py \
    /path_in/${INPUT} \
    /path_in/allTFs_hg38.txt \
    --method grnboost2 \
    --output /path_out/${OUTPUT} \
    --num_workers 130 \
    --seed 777