#include "./matrix_mul.h"
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>

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

// fixed: takes block_dim instead of sqrt(size)
void matrix_mul_add(double* A, double* B, double* C, int block_dim)
{
   for(int i=0; i < block_dim; i++)
   {
        for(int a=0; a < block_dim; a++)
        {
            double sum = 0;
            for(int b=0; b < block_dim; b++)
            {
                sum += (A[(i * block_dim) + b] * B[(b * block_dim) + a]);
            }
            C[(i * block_dim) + a] += sum;
        }
   }
}

// shift left
void shift_left(double *A, double *recv, int block_size, int i, int j, int q)
{
    int left  = i * q + ((j - 1 + q) % q);
    int right = i * q + ((j + 1) % q);

    MPI_Sendrecv(A, block_size, MPI_DOUBLE,
                 left, 0,
                 recv, block_size, MPI_DOUBLE,
                 right, 0,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

    for (int k = 0; k < block_size; k++)
        A[k] = recv[k];
}

// shift up
void shift_up(double *B, double *recv, int block_size, int i, int j, int q)
{
    int above = ((i - 1 + q ) % q) * q + j;
    int below = ((i + 1) % q) * q + j;

    MPI_Sendrecv(B, block_size, MPI_DOUBLE,
                 above, 0,
                 recv, block_size, MPI_DOUBLE,
                 below, 0,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

    for(int k = 0; k < block_size; k++)
        B[k] = recv[k];
}

int main(int argc, char *argv[])
{
    int id, p, len, q;
    char hostname[300];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    q = (int)sqrt(p);

    if(q*q != p) {
        if(id == 0) printf("ERROR: p must be a perfect square.\n");
        MPI_Finalize();
        return 1;
    }

    int block_dim  = N / q;          // ✓ safe
    int block_size = block_dim * block_dim;

    double *local_A = malloc(block_size * sizeof(double));
    double *local_B = malloc(block_size * sizeof(double));
    double *local_C = malloc(block_size * sizeof(double));

    for(int i = 0; i < block_size; i++)
        local_C[i] = 0.0;

    double *collect_C = NULL;
    double *A = NULL;
    double *B = NULL;

    MPI_Get_processor_name(hostname, &len);
    printf("Hostname: %s, rank: %d, size: %d\n", hostname, id, p);
    fflush(stdout);

    int row = id / q;
    int col = id % q;

    if(id == 0)
    {
        A = malloc(N*N*sizeof(double));
        B = malloc(N*N*sizeof(double));

        for(int i = 0; i < N*N; i++) {
            A[i] = rand() % MAX_RAND;
            B[i] = rand() % MAX_RAND;
        }

        double *send_A = malloc(p * block_size * sizeof(double));
        double *send_B = malloc(p * block_size * sizeof(double));

        for(int r = 0; r < q; r++)
        {
            for(int c = 0; c < q; c++)
            {
                int rank = r*q + c;
                int offset = rank * block_size;

                for(int i = 0; i < block_dim; i++)
                    for(int j = 0; j < block_dim; j++)
                    {
                        send_A[offset + i*block_dim + j] =
                            A[(r*block_dim + i)*N + (c*block_dim + j)];

                        send_B[offset + i*block_dim + j] =
                            B[(r*block_dim + i)*N + (c*block_dim + j)];
                    }
            }
        }

        MPI_Scatter(send_A, block_size, MPI_DOUBLE,
                    local_A, block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Scatter(send_B, block_size, MPI_DOUBLE,
                    local_B, block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        free(send_A);
        free(send_B);
    }
    else {
        MPI_Scatter(NULL, block_size, MPI_DOUBLE,
                    local_A, block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Scatter(NULL, block_size, MPI_DOUBLE,
                    local_B, block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    double *recv = malloc(block_size * sizeof(double));

    // Initial Cannon shifts
    for(int r = 0; r < row; r++)
        shift_left(local_A, recv, block_size, row, col, q);

    for(int c = 0; c < col; c++)
        shift_up(local_B, recv, block_size, row, col, q);

    double start = MPI_Wtime();

    for(int k = 0; k < q; k++)
    {
        matrix_mul_add(local_A, local_B, local_C, block_dim);

        shift_left(local_A, recv, block_size, row, col, q);
        shift_up(local_B, recv, block_size, row, col, q);
    }

    double end = MPI_Wtime();

    free(recv);

    if(id == 0)
        collect_C = malloc(N*N*sizeof(double));

    MPI_Gather(local_C, block_size, MPI_DOUBLE,
               collect_C, block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(local_A);
    free(local_B);
    free(local_C);

    if(id == 0)
    {
        double checksum = 0.0;
        for(int i = 0; i < N*N; i++)
            checksum += collect_C[i];

        printf("parallel check sum: %f\n", checksum);
        printf("Time Elapsed: %f seconds\n", end-start);
        fflush(stdout);

        free(collect_C);
        free(A);
        free(B);
    }

    MPI_Finalize();
    return 0;
}
