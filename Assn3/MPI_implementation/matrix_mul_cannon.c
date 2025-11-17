#include "./matrix_mul.h"
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>

//peusdo code for matrix multiplication found here:
//https://en.wikipedia.org/wiki/Matrix_multiplication_algorithm

//Input: matrices A and B
    //Let C be a new matrix of the appropriate size
    //For i from 1 to n:
        //For j from 1 to p:
            //Let sum = 0
            //For k from 1 to m:
                //Set sum ← sum + Aik × Bkj
            //Set Cij ← sum
    //Return C

//implementation of the puesdocode above
double** matrix_mul(double** A, double** B)
{
   double ** C = malloc(N*sizeof(double *)); 
   for(int i=0; i< N; i++)
   {
        C[i] = malloc(N*sizeof(double));
   }

   for(int i=0; i < N; i++)
   {
        for(int a=0; a < M; a++)
        {
            double sum = 0;
            for(int b=0; b < N; b++)
            {
                sum += (A[i][b] * B[b][a]);
            }
            C[i][a] = sum;
        }
   }

   return C;
}

//required for the cannons implementation and it now handles flattened data
void matrix_mul_add(double* A, double* B, double* C, int size)
{
   int n = (int)sqrt(size);
   for(int i=0; i < n; i++)
   {
        for(int a=0; a < n; a++)
        {
            double sum = 0;
            for(int b=0; b < n; b++)
            {
                sum += (A[(i * n) + b] * B[(b * n) + a]);
            }
            C[(i * n) + a] += sum;
        }
   }
}

//this will need to be converted into a usable 2d array however, but not just yet
//shift left, needs a recv buffer to handle the send and recv
void shift_left(double *A,  double *recv, int block_size, int i, int j, int q)
{
    //offset by row amount each process, and then move over to the column based on the which process is being used
    /*
    |
    |---> right and left are just columsn right next to the moving column
    */
    int left = i * q + ((j - 1 + q) % q);
    int right = i * q + ((j + 1) % q);  

    //covered briefly in the slides and would be more annoying to implement with just sends and recvs
    MPI_Sendrecv(A, block_size, MPI_DOUBLE, left, 0, recv, block_size, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    //can't have the same buffer so we have a new one to handle the recv
    for (int i = 0; i < block_size; i++)
    {
      A[i] = recv[i];  
    }
}

//shift up same as shift left but with rows now
void shift_up(double *B, double *recv, int block_size, int i, int j, int q)
{
    //ofset by processes because bottom row has to be sent to top
    int above = ((i - 1 + q ) % q) * q + j;

    //this calculation is easier because it just drops one lower
    int below = ((i + 1) % q) * q + j;

    MPI_Sendrecv(B, block_size, MPI_DOUBLE, above, 0, recv, block_size, MPI_DOUBLE, below, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for(int i = 0; i < block_size; i++)
    {
        B[i] = recv[i];
    }

}

int main(int argc, char *argv[])
{
    //MPI variables
    int id, p, len, q, rows, cols;
    char hostname[300];
    
    //set seed so runs are predictable
    srand(5050);

    //standard setup for MPI programs
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    //https://docs.open-mpi.org/en/v5.0.x/man-openmpi/man3/MPI_Get_processor_name.3.html#mpi-get-processor-name
    MPI_Get_processor_name(hostname, &len);
    printf("Hostname: %s, rank: %d, size: %d\n", hostname, id, p);
    fflush(stdout);

    q=(int)sqrt(p);

    //local arrays for each process
    double * local_A = malloc((N/q)*(N/q)*sizeof(double));
    double * local_B = malloc((N/q)*(N/q)*sizeof(double));
    double * local_C = malloc((N/q)*(N/q)*sizeof(double));
    for (int i = 0; i < (N/q)*(N/q); i++)
    {
        local_C[i] = 0.0;
    }

    //collect array for root
    double *collect_C;
    double * A;
    double * B;
    //setting up i, j so we can do the other operations later

    rows = id / q;
    cols = id % q;

    //what the root allocates blocks via connon's method, handling 2d blocks, but we need to work with them flattened first
    //https://en.wikipedia.org/wiki/Cannon%27s_algorithm
    if(id == 0)
    {
        A = malloc(N*N*sizeof(double));
        B = malloc(N*N*sizeof(double));
        //setup the matrixes 
        for(int i = 0; i < N; i++)
        {
            for(int a = 0; a < N; a++)
            {
                A[(i * N) + a] = (double)(rand() % MAX_RAND);
                B[(i * N) + a] = (double)(rand() % MAX_RAND);
            }
        }

        double *send_A = malloc(p * (N/q) * (N/q) * sizeof(double));
        double *send_B = malloc(p * (N/q) * (N/q) * sizeof(double));


        //give out sections of the matrices
        for(int i = 0; i < q; i++)
        {
            for(int a = 0; a < q; a++)
            {
                int rank = i * q + a;
                int element = rank * (N/q) * (N/q);

                int block_size = (N/q) * (N/q);

                for(int b = 0; b < N/q; b++)
                {
                    for(int d = 0; d < N/q; d++)
                    {
                        send_A[element + b *(N/q) + d] = A[(i*(N/q)+b)*N + (a*(N/q)+d)];
                        send_B[element + b *(N/q) + d] = B[(i*(N/q)+b)*N + (a*(N/q)+d)];
                    }
                }
                
            }
        }

         //send to the specificed process
                MPI_Scatter(send_A, (N/q) * (N/q), MPI_DOUBLE, local_A, (N/q) * (N/q), MPI_DOUBLE, 0, MPI_COMM_WORLD);
                MPI_Scatter(send_B, (N/q) * (N/q), MPI_DOUBLE, local_B, (N/q) * (N/q), MPI_DOUBLE, 0, MPI_COMM_WORLD);
            
                free(send_A);
                free(send_B);
    } else {
        MPI_Scatter(NULL, (N/q)*(N/q), MPI_DOUBLE, local_A, (N/q)*(N/q), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(NULL, (N/q)*(N/q), MPI_DOUBLE, local_B, (N/q)*(N/q), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    //Shift A(i,j) left by i positions (mod q)
    double *recv = malloc((N/q) * (N/q) * sizeof(double));
    for(int r = 0; r < rows; r++)
    {
        shift_left(local_A,recv, (N/q)*(N/q), rows, cols, q);
    }

    //Shift B(i,j) up by j positions (mod q)
    for(int c = 0; c < cols; c++)
    {
        shift_up(local_B, recv, (N/q)*(N/q), rows, cols, q);
    }

    // Step 2: Compute
    //Cij = 0 already handled by the zeroed out 2d C array

    //the meat and potatoes of the program/the parallel section we need to time
    double start = MPI_Wtime();
    for(int k = 0; k < q; k++)
    {
        matrix_mul_add(local_A, local_B, local_C, (N/q)*(N/q));   // Local block multiplication how do we accumulate this

        shift_left(local_A, recv, (N/q)*(N/q), rows, cols, q);
        shift_up(local_B, recv, (N/q)*(N/q), rows, cols, q);
    }
    double end = MPI_Wtime();

    free(recv);

    if(id == 0)
    {
        //initialize the C_local array
        collect_C = malloc(N*N*sizeof(double));
    }
    // Step 3: Gather results if needed Assemble C from all Cij blocks
    MPI_Gather(local_C, (N/q)*(N/q), MPI_DOUBLE, collect_C, (N/q)*(N/q), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);

    free(local_C);
    free(local_A);
    free(local_B);

    if(id == 0)
    {
        //For small values of N create a sequential run and then compare checksum

        double checksum_para = 0.0;
        for (int i = 0; i < N*N; i++) {
            checksum_para += collect_C[i];
        }

        free(collect_C);

        if(N < 200)
        {
            //convert to 2d array to run the check on it
            double **A2 = malloc(N * sizeof(double*));
            double **B2 = malloc(N * sizeof(double*));
            for (int i = 0; i < N; i++)
            {
                A2[i] = malloc(N * sizeof(double));
                B2[i] = malloc(N * sizeof(double));
            }
            for(int a = 0; a < N; a++)
            {
              for (int j = 0; j < N; j++)
                {
                    A2[a][j] = A[a*N + j];
                    B2[a][j] = B[a*N + j];
                }
            }

            double** C_check = matrix_mul(A2, B2);
        
            double checksum_seq = 0.0;

            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    checksum_seq += C_check[i][j];
                }
            }

            if(checksum_seq == checksum_para)
            {
                printf("check sum passed same as sequential version\n");
            }
            
            for (int i = 0; i < N; i++) {
                free(A2[i]);
                free(B2[i]);
            }
            free(A2);
            free(B2);
            free(A);
            free(B);

        }

        printf("parallel check sum: %f\n", checksum_para);
        printf("Time Elasped: %f seconds\n", end-start);
        
    }
    MPI_Finalize();
    return 0;
}