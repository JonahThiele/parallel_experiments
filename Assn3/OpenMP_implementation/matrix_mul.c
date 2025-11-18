#include "./matrix_mul.h"
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#define MAX_RAND 100  

//standard matrix muliplication
double** matrix_mul(double** A, double** B, int N)
{
    double ** C = malloc(N * sizeof(double *)); 
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
    if(argc < 3)
    {
        printf("add N and number of threads");
        return -1;
    }

    //get N and number of threads
    int N = atoi(argv[1]);
    int n_threads = atoi(argv[2]);

    omp_set_num_threads(n_threads);

    //we need to print the number of threads
    #pragma omp parallel
    {
        #pragma omp single
        printf("%d threads\n", omp_get_num_threads());
    }

    // create the matrixes
    double **A = malloc(N*sizeof(double*));
    double **B = malloc(N*sizeof(double*));
    double **C = malloc(N*sizeof(double*));
    for(int i = 0; i < N; i++)
    {
        A[i] = malloc(N*sizeof(double));
        B[i] = malloc(N*sizeof(double));
        C[i] = malloc(N*sizeof(double));
    }

    //load the matrices
    srand(5050);
    for(int i = 0; i < N; i++)
    {
        for(int a = 0; a < N; a++)
        {
            A[i][a] = rand() % MAX_RAND;
            B[i][a] = rand() % MAX_RAND;
            C[i][a] = 0;
        }
    }

    double start = omp_get_wtime();
    #pragma omp parallel for schedule(dynamic, 4)
    //#pragma omp parallel for schedule(static, 10)
    for(int i = 0; i < N; i++)
    {
        for(int a = 0; a < N; a++)
        {
            double sum = 0;
            for(int j = 0; j < N; j++)
            {
                sum += A[i][j] * B[j][a];
            }
            C[i][a] = sum;
        }
    }
    double end = omp_get_wtime();

    // Sequential check for correctness
     double checksum_para = 0;
    for(int i = 0; i < N; i++)
    {
        for(int a = 0; a < N; a++)
        {
            checksum_para  += C[i][a];
        }
    }
    if(N < 200)
    {
        double **C_check = matrix_mul(A, B, N);
        double checksum_seq  = 0;
        for(int i = 0; i < N; i++)
        {
            for(int a = 0; a < N; a++)
            {
                checksum_seq  += C_check[i][a];
            }
        }

        if(checksum_para != checksum_seq)
        {
            printf("Check sum failed.\n");
        }else{
            printf("Check sum succeed\n");
        }
        for(int i = 0; i < N; i++)
        {
            free(C_check[i]);
        }
        free(C_check);
    }
    printf("checksum: %f\n", checksum_para);
    printf("total time: %f seconds\n", end - start);

    // clean up the mess
    for(int i = 0; i < N; i++)
    {
        free(A[i]);
        free(B[i]);
        free(C[i]);
    }
    free(A); 
    free(B); 
    free(C); 

    return 0;
}
