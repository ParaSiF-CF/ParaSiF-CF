#!/bin/bash

# Slurm job options (job-name, compute nodes, job time)
#SBATCH --job-name=ParaSiF_CF_t01
#SBATCH --time=10:20:00
#SBATCH --exclusive
#SBATCH --export=none
#SBATCH --account=***

#SBATCH --partition=standard
#SBATCH --qos=standard
# #SBATCH --reservation=shortqos

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1

#SBATCH hetjob

#SBATCH --partition=standard
#SBATCH --qos=standard
# #SBATCH --reservation=shortqos

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

pwd; hostname; date

# Setup the job environment (this module needs to be loaded before any other modules)
module load epcc-job-env

module -s restore /etc/cray-pe.d/PrgEnv-gnu

. config/ParaSiF_CF_gcc.conf

module list

pip list

# ./SCRIPTS/compileCS

domainFluid=${PWD}/RESU/fluidDomain
domainStructure=${PWD}/RESU/structureDomain

# Get application name for fluid and structure domain
solverFluid=./cs_solver
solverStructure=structureDomainRun.py

export OMP_NUM_THREADS=1

srun  --distribution=block:block --hint=nomultithread --het-group=1 --cpu-bind=rank --chdir=${domainFluid} ${solverFluid} --mpi :\
     --distribution=block:block --hint=nomultithread --het-group=0 --cpu-bind=rank --chdir=${domainStructure} python3 ${solverStructure}

