
#!/bin/bash
#PBS -l select=1:ncpus=16:mem=96gb:ngpus=4:gpu_type=RTX6000
#PBS -l walltime=24:00:00

tid=$PBS_ARRAY_INDEX
module load anaconda3/personal
source activate snRNAseq

cellbender remove-background \
  --input /rds/general/project/lms-scott-raw/live/combined_ngs_cellranger_outs_human_WAT_obese_lean_wl/Processed_combined_introns/Pool1OB/outs/raw_feature_bc_matrix.h5 \
  --output /rds/general/project/lms-scott-raw/live/snRNA_seq_processing/human_WAT/Phase_4/cellbender/cellbender.Pool1B.output.h5 \
  --cuda \
  --expected-cells 5000 \
  --total-droplets-included 30000 \
  --fpr 0.01 \
  --epochs 150

cellbender remove-background \
  --input /rds/general/project/lms-scott-raw/live/combined_ngs_cellranger_outs_human_WAT_obese_lean_wl/Processed_combined_introns/WS4/Pool1WL/raw_feature_bc_matrix.h5 \
  --output /rds/general/project/lms-scott-raw/live/snRNA_seq_processing/human_WAT/Phase_4/cellbender/cellbender.WS4.output.h5 \
  --cuda \
  --expected-cells 5000 \
  --total-droplets-included 30000 \
  --fpr 0.01 \
  --epochs 150

cellbender remove-background \
  --input /rds/general/project/lms-scott-raw/live/combined_ngs_cellranger_outs_human_WAT_obese_lean_wl/Processed_combined_introns/WS5/Pool1LN/raw_feature_bc_matrix.h5 \
  --output /rds/general/project/lms-scott-raw/live/snRNA_seq_processing/human_WAT/Phase_4/cellbender/cellbender.WS5.output.h5 \
  --cuda \
  --expected-cells 5000 \
  --total-droplets-included 30000 \
  --fpr 0.01 \
  --epochs 150
