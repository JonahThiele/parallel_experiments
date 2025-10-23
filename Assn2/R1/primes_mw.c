#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define WORK 100
#define RESULT 101
#define STOP 102

//simple naive trial divison to find primes
int n_primes(int start, int end)
{
    int primes = 0;
    //check the entire range given
    for(int i = start; i <= end; i++)
    {
        bool prime = true;
        //prime if 2
        if(i == 2)
        {
            primes++;
            continue;
        }
        //don't check evens
        else if( i % 2 == 0)
        {
            prime = false;
        }
        else{
            //check divisors up to square root
            for(int a = 3; a * a <= i; a += 2)
            {
                if(i % a == 0)
                {
                    prime = false;
                    break;
                }
            }
            if(prime)
            {
                primes++;
            }
        }
    }
    return primes;
}

int main (int argc, char *argv[])
{

    int i;

    //comm variables
    int id;
    int p;

    //arguments 
    int chunk;
    int n;

    //workers
    int primes = 0;
    int work[2];

    //root
    int total_primes = 0;

    //show hostname
    char hostname[300];
    int len;

    //timing 
    double start, end;

    //strtol code reference from here: https://stackoverflow.com/questions/9748393/how-can-i-get-argv-as-int
    if(argc < 3)
    {
        printf("the two required arguments N and CHUNK where not supplied\n");
        fflush(stdout);
        return -1;
    }

    //handle n and chunk args
    n = strtol(argv[1], NULL, 10);
    chunk = strtol(argv[2], NULL, 10);

    //make sure chunk is not 0
    if(chunk < 1)
    {
        printf("chunk size is too small must be greater than 0");
        return -1;
    }

    //make sure n is not the same as 2
    if(n < 2)
    {
        printf("n is too small it must be greater than 2");
        return -1;
    }

    //MPI standard setup stuff
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    //using this instead because its more portable https://stackoverflow.com/questions/22174901/gethostname-function-missing-in-openmpi
    //https://docs.open-mpi.org/en/v5.0.x/man-openmpi/man3/MPI_Get_processor_name.3.html#mpi-get-processor-name
    MPI_Get_processor_name(hostname, &len);
    printf("Hostname: %s, rank: %d, size: %d\n", hostname, id, p);
    fflush(stdout);

    //if master start handing out work
    if (id == 0)
    {
        start = MPI_Wtime();
        printf("N=%d, P=%d, CHUNK=%d\n", n, p, chunk);
        //cylic allocation but allows for lots of idle time
        for(int i = 0; i * chunk < n; i++)
        {
            int start = i * chunk + 2;
            int end = start + chunk - 1;
            if(start > n) {
                break;  
            }
            if(end > n)
            {
                end = n;
            }
            printf("master send [%d, %d] to worker %d hostname:\n", start, end, (i % (p - 1) + 1), hostname);
            work[0] = start;
            work[1] = end;
            MPI_Send(work, 2, MPI_INT, (i % (p - 1)) + 1, WORK, MPI_COMM_WORLD);
            //idle time is here waiting for reponse from worker
            MPI_Recv(&primes, 1, MPI_INT, (i % (p - 1)) + 1, RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_primes += primes; 
        }

        //set out the stop to all workers
        for(int i = 1; i < p; i++)
        {
            MPI_Send(work, 2, MPI_INT, i, STOP, MPI_COMM_WORLD);
        }
        end = MPI_Wtime();
    } else {
        //recv work and send back results until recv stop
        while(true)
        {
            MPI_Status status;
            MPI_Recv(work, 2, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if(status.MPI_TAG == STOP)
            {
                break;
            }
            primes = n_primes(work[0], work[1]);
            printf("worker %d: processed [%d..%d] -> %d primes\n",id, work[0], work[1], primes);
            fflush(stdout);
            MPI_Send(&primes, 1, MPI_INT, 0, RESULT, MPI_COMM_WORLD);
        }
    }
    //only check when n is really small
    if(id == 0 && n <= 200000)
    {
        //sequential check same code as above
        int check_primes = 0;
        for(int i = 2; i <= n; i++)
        {
            bool prime = true;
            if ( i == 2)
            {
                check_primes++;
                continue;
            }
            else if( i % 2 == 0)
            {
                prime = false;
            }
            else{
                for(int a = 3; a * a <= i; a += 2)
                {
                    if(i % a == 0)
                    {
                        prime = false;
                        break;
                    }
                }
                if(prime)
                {
                    check_primes++;
                }
            }
        }
        //handle a failed check
        if(check_primes != total_primes)
        {
            printf("failed the correctness check -> ERROR <-\n workers: %d seq: %d\n", total_primes, check_primes);
            fflush(stdout); 
        }
    }
    //final summary
    if(id == 0)
    {
        printf("Final summary from master:\n");
        printf("N=%d, P=%d, CHUNK=%d\n", n, p, chunk);
        printf("total primes = %d\n", total_primes);
        printf("elapsed time = %fs\n", end-start);
        fflush(stdout);
    }
    MPI_Finalize();
    return 0;
}