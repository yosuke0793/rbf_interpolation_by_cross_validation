#!/usr/bin/env sh
#----------------------------------------------#
#       TORQUE BATCH JOB SCRIPT TEMPLATE       #
#                     v6.0                     #
#----------------------------------------------#
#PBS -N rbfinter
#PBS -l nodes=1:ppn=36:iruka
#PBS -l walltime=24:00:00
#PBS -e ./stderr/
#PBS -o ./stdout/
#----------------------------------------------#
cd $PBS_O_WORKDIR
rm -rf std*
mkdir -p stdout stderr
[ $PBS_O_HOST == eggplant ] && opt="-bind-to-none"
[ $PBS_O_HOST == orca     ] && opt="-bind-to none"
nprocs=`wc -l < $PBS_NODEFILE`
export OMP_NUM_THREADS=2
#----------------------------------------------#
rm -f ./inputs/dataset.txt
rm -f ./results/*
cp  -f ../datamake/results/dataset.txt ./inputs/
mpirun -machinefile $PBS_NODEFILE -np 18 $opt rbf_interpolation 
python ./visualize_results.py
