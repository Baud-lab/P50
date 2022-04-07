#!/bin/bash
#$ -N map
#$ -o /users/abaud/abaud/micturition/output/output/$JOB_NAME_$TASK_ID.out
#$ -e /users/abaud/abaud/micturition/output/error/$JOB_NAME_$TASK_ID.err
#$ -t 1
#$ -q long-sl7

# -l virtual_free=50G

python /users/abaud/abaud/micturition/code/GWAS/map_DGE_LOCO_GPKronSum.py $SGE_TASK_ID exp_DGA_7_Dec_2021 citg

#submit with qsub < /users/abaud/abaud/micturition/code/GWAS/qsub_map_DGE_LOCO.sh
