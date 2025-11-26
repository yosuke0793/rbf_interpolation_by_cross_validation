rm -f ./inputs/dataset.txt
cp  -f ../datamake/results/dataset.txt ./inputs/
### Single processor run
# rbf_interpolation
### Parallel run
# nprocs = 4
mpirun -machinefile $PBS_NODEFILE -np $nprocs $opt rbf_interpolation 
python ./visualize_results.py
