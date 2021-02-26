#PBS -l nodes=1:ppn=28
#PBS -l walltime=48:00:00
#PBS -S /bin/bash
#PBS -N wsgups
#PBS -j oe
#PBS -o LOG

# this bash script runs the pipeline in an HPC environment
# please adapt to your particular HPC settings

cd $PBS_O_WORKDIR

source activate base

snakemake --cluster "qsub -q short" --use-conda --local-cores 28 -j 28 -k
