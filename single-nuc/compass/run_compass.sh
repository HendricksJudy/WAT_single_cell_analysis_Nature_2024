#!/usr/bin/bash

#PBS -lselect=1:ncpus=80:mem=600gb
#PBS -lwalltime=48:00:00

export PATH=$Home/rds/general/user/username/home/.conda/envs/personal/bin/python:$PATH
cd /compass/log1p_mean_adipocytes # change accordingly

module load anaconda3/personal
unset PYTHONPATH
source activate py_compass

compass --data expression_matrix_log1p_per_sample_mean.tsv --num-processes 50 --species homo_sapiens