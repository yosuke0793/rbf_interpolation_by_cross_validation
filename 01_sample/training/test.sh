#!/usr/bin/env sh
#----------------------------------------------#
#       TORQUE BATCH JOB SCRIPT TEMPLATE       #
#                     v6.0                     #
#----------------------------------------------#
#PBS -N job_name
#PBS -l nodes=1:ppn=1:tan
#PBS -l walltime=24:00:00
#PBS -e ./stderr/
#PBS -o ./stdout/
#----------------------------------------------#
cd $PBS_O_WORKDIR
mkdir -p stdout stderr
[ $PBS_O_HOST == eggplant ] && opt="-bind-to-none"
[ $PBS_O_HOST == orca     ] && opt="-bind-to none"
nprocs=`wc -l < $PBS_NODEFILE`
export OMP_NUM_THREADS=$nprocs
#----------------------------------------------#

./a.out
#mpirun -machinefile $PBS_NODEFILE -np $nprocs $opt ./a.out

