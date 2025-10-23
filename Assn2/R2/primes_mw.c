#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

//setup the tags to use
#define WORK 100
#define RESULT 101
#define STOP 102
#define REQUEST 103

//prime code from R1
int n_primes(int start, int end)
{
    int primes = 0;
    for(int i = start; i <= end; i++)
    {
        bool prime = true;
        if(i == 2)
        {
            primes++;
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
                primes++;
            }
        }
    }
    return primes;
}

int main (int argc, char *argv[])
{

    //global copies for MPI info
    int id;
    int p;
    int chunk;
    int n;
    
    //for workers
    int primes = 0;
    int total_primes = 0;

    //show hostname
    char hostname[300];
    int len;

    //worker range
    int work[2];

    //timing 
    double start, end;

    //strtol code reference from here: https://stackoverflow.com/questions/9748393/how-can-i-get-argv-as-int
    if(argc < 3)
    {
        printf("the two required arguments N and CHUNK where not supplied\n");
        return -1;
    }
    //get n and chunk args
    n = strtol(argv[1], NULL, 10);
    chunk = strtol(argv[2], NULL, 10);

    //chunk can't be zero
    if(chunk < 1)
    {
        printf("chunk size is too small must be greater than 0\n");
        return -1;
    }

    //n can't be 2
    if(n < 2)
    {
        printf("n is too small it must be greater than 2\n");
        return -1;
    }

    //standard MPI setup
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    //using this instead because its more portable https://stackoverflow.com/questions/22174901/gethostname-function-missing-in-openmpi
    //https://docs.open-mpi.org/en/v5.0.x/man-openmpi/man3/MPI_Get_processor_name.3.html#mpi-get-processor-name
    MPI_Get_processor_name(hostname, &len);

    printf("Hostname: %s, rank: %d, size: %d\n", hostname, id, p);
    fflush(stdout);

    //use dynamic allocation of work
    if (id == 0)
    {
        start = MPI_Wtime();
        printf("N=%d, P=%d, CHUNK=%d\n", n, p, chunk);
        fflush(stdout);
        int i = 0;
        int workers_working = p - 1;
        //see if workers are still processing ranges
        while(workers_working > 0)
        {
            MPI_Status status;
            MPI_Recv(&primes, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if(status.MPI_TAG == REQUEST)
            {
                //we still have ranges to work on
                if(i * chunk < n)
                {
                    int start = i * chunk + 2;
                    int end = start + chunk - 1; //no overlap
                    //make sure start doesn't come before end
                    if(start > n) {
                        break;  
                    } 
                    //handle the last chunk
                    if(end > n)
                    {
                        end = n;
                    }
                    printf("master send [%d, %d] to worker %d\n", start, end, (i % (p - 1) + 1));
                    fflush(stdout);
                    work[0] = start;
                    work[1] = end;
                    MPI_Send(work, 2, MPI_INT, status.MPI_SOURCE, WORK, MPI_COMM_WORLD);
                    i++;
                } else {
                    //send  STOP because of no more work
                    MPI_Send(work, 2, MPI_INT, status.MPI_SOURCE, STOP, MPI_COMM_WORLD);
                    workers_working--;
                }
            }else{
                //sum up the results
                total_primes += primes;
            }
        }
        end = MPI_Wtime();
    } else {
        while(true)
        {
            //first send request
            MPI_Send(&primes, 1, MPI_INT, 0, REQUEST, MPI_COMM_WORLD);

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
        //sequential check same as R1
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
        //handle the works and sequential check not matching
        if (check_primes != total_primes)
        {
            printf("failed the correctness check -> ERROR <-\n workers: %d seq: %d\n", total_primes, check_primes);
            fflush(stdout);
        }
    }
    //summary
    if(id == 0){
        printf("Final summary from master:\n");
        printf("N=%d, P=%d, CHUNK=%d\n", n, p, chunk);
        printf("total primes = %d\n", total_primes);
        printf("elapsed time = %fs\n", end-start);
        fflush(stdout);
    }
    MPI_Finalize();
    return 0;
}