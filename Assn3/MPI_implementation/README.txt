compile the code as such:

mpicc matrix_mul.c -o mul 

run the code with:

mpirun --machinefile machinefile -np {number of process} ./mul {size of N} {number of threads}