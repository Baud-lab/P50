#BSUB -o /homes/abaud/P50_HSrats/output/output/%J_%I.out
#BSUB -e /homes/abaud/P50_HSrats/output/error/%J_%I.err
#BSUB -M 10000
#BSUB -R "rusage[mem=20000]"
#BSUB -J "VD[1-789]"
#change 789 above to total number of columns in metabolite table in HDF5 file
#BSUB -q highpri

python /homes/abaud/P50_HSrats/code/variance_decomposition/metabo/exp_varianceDecomp_LOCO.py $LSB_JOBINDEX pruned_dosages

#submit with bsub < /homes/abaud/P50_HSrats/code/variance_decomposition/metabo/qsub_varianceDecomp_LOCO.sh
