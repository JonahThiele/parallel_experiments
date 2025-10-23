## how to build R1 (make sure you are in R1 directory)
mpicc primes_mw.c -o primes_mw

## how to build R2 (make sure you are in R2 directory)
mpicc primes_mw.c -o primes_mw

## how to run R1 without hostfile
mpirun --np (the number of processors) ./primes_mw N(the number of primes) CHUNK(the size of chunks)

## how to run R1 with hostfile
mpirun --machinefile N(name of machine file) --np (the number of processors) ./primes_mw N(the number of primes) CHUNK(the size of chunks)

## how to run R2 without hostfile
mpirun --np (the number of processors) ./primes_mw N(the number of primes) CHUNK(the size of chunks)

## how to run R2 with hostfile
mpirun --machinefile (name of machine file) --np (the number of processors) ./primes_mw N(the number of primes) CHUNK(the size of chunks)

## arguments
N -> the range of numbers from 2 to N that primes will be found in
CHUNK -> the size of the range to find primes that will be given to each process

## limits or assumptions
N and CHUNK are integers, that fit within the standard ranges of C integers
CHUNK is assumed smaller than N

## primality test method
the method used is just a naive trial division by 