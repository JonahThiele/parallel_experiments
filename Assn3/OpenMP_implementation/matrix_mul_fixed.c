#include "./matrix_mul.h"
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#define MAX_RAND 100  // example max random value

double** matrix_mul(double** A, double** B, int N)
{
    double ** C = malloc(N * sizeof(double *)); 
    for(int i = 0; i < N; i++)
        C[i] = malloc(N * sizeof(double));

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
    if(argc < 3)
    {
        printf("Usage: %s <matrix_size> <num_threads>\n", argv[0]);
        return -1;
    }

    int N = atoi(argv[1]);
    int n_threads = atoi(argv[2]);

    // Set number of threads
    omp_set_num_threads(n_threads);

    // Print number of threads inside a parallel region
    #pragma omp parallel
    {
        #pragma omp single
        printf("Running with %d threads\n", omp_get_num_threads());
    }

    // Allocate matrices
    double **A = malloc(N * sizeof(double *));
    double **B = malloc(N * sizeof(double *));
    double **C = malloc(N * sizeof(double *));
    for(int i = 0; i < N; i++)
    {
        A[i] = malloc(N * sizeof(double));
        B[i] = malloc(N * sizeof(double));
        C[i] = malloc(N * sizeof(double));
    }

    // Initialize matrices
    srand(5050);
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            A[i][j] = rand() % MAX_RAND;
            B[i][j] = rand() % MAX_RAND;
            C[i][j] = 0.0;
        }
    }

    // Parallel matrix multiplication
    double start = omp_get_wtime();
    #pragma omp parallel for schedule(dynamic, 4)
    for(int i = 0; i < N; i++)
    {
        for(int k = 0; k < N; k++)
        {
            double sum = 0.0;
            for(int j = 0; j < N; j++)
            {
                sum += A[i][j] * B[j][k];
            }
            C[i][k] = sum;
        }
    }
    double end = omp_get_wtime();

    // Sequential check for correctness
    double **C_check = matrix_mul(A, B, N);

    double checksum_para = 0.0;
    double checksum_seq  = 0.0;
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            checksum_para += C[i][j];
            checksum_seq  += C_check[i][j];
        }
    }

    if(checksum_para == checksum_seq)
        printf("Check sum passed! Parallel result matches sequential.\n");
    else
        printf("Check sum failed!\n");

    printf("Elapsed time: %f seconds\n", end - start);

    // Free memory
    for(int i = 0; i < N; i++)
    {
        free(A[i]);
        free(B[i]);
        free(C[i]);
        free(C_check[i]);
    }
    free(A); free(B); free(C); free(C_check);

    return 0;
}
