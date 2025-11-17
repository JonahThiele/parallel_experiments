#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define MAX_RAND 100

//use flat arrays because its a lot easier to send flat arrays
void matrix_mul_flat(double *A, double *B, double *C, int N, int rows)
{
    for(int i = 0; i < rows; i++)
    {
        for(int a = 0; a < N; a++)
        {
            double sum = 0;
            for(int b = 0; b < N; b++)
            {
                sum += A[i*N + b] * B[b*N + a];
            }

            C[i*N + a] = sum;
        }
    }
}

//standard matrix multiplication with 2d arrays
double** matrix_mul(double** A, double** B, int N)
{
    double **C = malloc(N * sizeof(double *));
    for(int i = 0; i < N; i++)
    {
        C[i] = malloc(N * sizeof(double));
    }

    for(int i = 0; i < N; i++)
    {
        for(int a = 0; a < N; a++)
        {
            double sum = 0;
            for(int b = 0; b < N; b++)
            {
                sum += A[i][b] * B[b][a];
            }
            C[i][a] = sum;
        }
    }

    return C;
}

int main(int argc, char* argv[])
{
    int id;
    int p;
    char hostname[200];
    int len;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if(argc < 2)
    {
        if(id == 0) printf("Add the values of N\n");
        MPI_Finalize();
        return -1;
    }

    MPI_Get_processor_name(hostname, &len);
    printf("Hostname: %s, rank: %d, size: %d\n", hostname, id, p);
    fflush(stdout);

    int N = atoi(argv[1]);
    int rows = N / p;
    int remainder = N % p;
    int local_rows = rows;

    if(id == 0 && remainder > 0)
    {
        local_rows += remainder;
    }

    double *local_A = malloc(local_rows * N * sizeof(double));
    double *local_C = malloc(local_rows * N * sizeof(double));
    double *B = malloc(N*N*sizeof(double));

    double *handle_remainder_A = NULL;
    double *handle_remainder_C = NULL;
    double *A = NULL;
    double *C = NULL;

    if(id == 0)
    {
        A = malloc(N*N*sizeof(double));
        C = malloc(N*N*sizeof(double));
        handle_remainder_A = malloc((N - remainder)*N*sizeof(double));
        handle_remainder_C = malloc((N - remainder)*N*sizeof(double));
        srand(5050);
        for(int i = 0; i < N*N; i++)
        {
            A[i] = rand() % MAX_RAND;
            B[i] = rand() % MAX_RAND;
            C[i] = 0;
        }

        for(int i = 0; i < N - remainder; i++)
        {
            for(int a = 0; a < N; a++)
            {
                handle_remainder_A[i * N + a] = A[(i + remainder) * N + a];
                handle_remainder_C[i * N + a] = 0;
            }
        }
    }

    // give B to all processes
    MPI_Bcast(B, N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // distribute rows using MPI_Send / MPI_Recv
    if(id == 0)
    {
        // root keeps its own local_A
        for(int i = 0; i < local_rows; i++)
            for(int j = 0; j < N; j++)
                local_A[i*N + j] = handle_remainder_A[i*N + j];

        // send chunks to other processes
        for(int i = 1; i < p; i++)
        {
            MPI_Send(handle_remainder_A + i*rows*N, rows*N, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        MPI_Recv(local_A, rows*N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    double start = MPI_Wtime();
    matrix_mul_flat(local_A, B, local_C, N, local_rows);
    double end = MPI_Wtime();

    // handle remainder rows on root
    if (id == 0 && remainder > 0)
    {
        matrix_mul_flat(A, B, C, N, remainder);
    }

    // gather results using MPI_Send / MPI_Recv
    if(id == 0)
    {
        // copy root's results
        for(int i = 0; i < local_rows; i++)
            for(int j = 0; j < N; j++)
                handle_remainder_C[i*N + j] = local_C[i*N + j];

        // receive from other processes
        for(int i = 1; i < p; i++)
        {
            MPI_Recv(handle_remainder_C + i*rows*N, rows*N, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    else
    {
        MPI_Send(local_C, rows*N, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }

    // combine results
    if(id == 0)
    {
        for(int i = 0; i < N - remainder; i++)
            for(int a = 0; a < N; a++)
                C[(i + remainder)*N + a] = handle_remainder_C[i*N + a];

        double checksum_para = 0;
        for(int i = 0; i < N*N; i++) checksum_para += C[i];

        if(N < 200)
        {
            double **A_test = malloc(N*sizeof(double*));
            double **B_test = malloc(N*sizeof(double*));
            for(int i = 0; i < N; i++)
            {
                A_test[i] = malloc(N*sizeof(double));
                B_test[i] = malloc(N*sizeof(double));
            }

            for(int i = 0; i < N; i++)
                for(int a = 0; a < N; a++)
                {
                    A_test[i][a] = A[i*N + a];
                    B_test[i][a] = B[i*N + a];
                }

            double **C_test = matrix_mul(A_test, B_test, N);
            double checksum_seq = 0;
            for(int i = 0; i < N; i++)
                for(int a = 0; a < N; a++)
                    checksum_seq += C_test[i][a];

            if(checksum_para != checksum_seq)
                printf("Checksum failed! the two matrices are not equal\n");

            printf("Checksum sequential: %f\n", checksum_seq);
        }

        printf("Checksum total: %f\n", checksum_para);
        printf("Parallel time: %lf secs\n", end-start);

        free(A);
        free(C);
        free(handle_remainder_A);
        free(handle_remainder_C);
    }

    free(B);
    free(local_A);
    free(local_C);

    MPI_Finalize();
    return 0;
}
