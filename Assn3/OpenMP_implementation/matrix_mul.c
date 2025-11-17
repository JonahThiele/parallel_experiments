#include "./matrix_mul.h"
#include "stdlib.h"
#include "stdio.h"
#include "omp.h"

double** matrix_mul(double** A, double** B, int N)
{
   double ** C = malloc(N*sizeof(double *)); 
   for(int i=0; i< N; i++)
   {
        C[i] = malloc(N*sizeof(double));
   }

   for(int i=0; i < N; i++)
   {
        for(int a=0; a < N; a++)
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

int main(int argc, char* argv[])
{   
    if (argc < 2)
    {
        printf("didn't enter the size fo the array\n");
        return -1;
    }else if (argc < 3)
    {
        printf("didn't enter the number of threads\n");
        return -1;
    }

    int N = atoi(argv[1]);
    int n_threads = atoi(argv[2]);

    //set the number of threads from the commandline
    omp_set_num_threads(n_threads);
    printf("Running with %d threads\n", omp_get_num_threads());
    double ** A = malloc(N*sizeof(double *));
    double ** B = malloc(N*sizeof(double *));
    double ** C = malloc(N*sizeof(double *)); 
    for(int i=0; i< N; i++)
    {
        C[i] = malloc(N*sizeof(double));
        B[i] = malloc(N*sizeof(double));
        A[i] = malloc(N*sizeof(double));
    }

    srand(5050);

    for(int i =0;i < N ;i++)
    {
        for(int a=0; a < N;a++)
        {
            A[i][a] = rand()%MAX_RAND;
            B[i][a] = rand()%MAX_RAND;
            C[i][a] = 0;
        }
    }

    double start = omp_get_wtime();

    //#pragma omp parallel for schedule(static) different scheduling strategy
    #pragma omp parallel for shared(A, B, C) schedule(dynamic, 4)
    for(int i = 0; i < N; i++)
    {   
        for(int k = 0; k < N; k++)
        {
            double sum = 0;
            for(int j = 0; j < N; j++)
            {
                sum += (A[i][j] * B[j][k]);
            }
            C[i][k] = sum;
        }
    }
    double end = omp_get_wtime();


    //check sum stuff

    double** C_check = matrix_mul(A, B, N);    
    double checksum_seq = 0.0;
    double checksum_para = 0.0;

    for(int i = 0; i < N; i++)
    {
        for(int a = 0; a < N; a++)
        {
            checksum_para += C[i][a];
            checksum_seq += C_check[i][a];
        }
    }

    if(checksum_seq == checksum_para)
    {
        printf("check sum passed same as sequential version\n");
    }
            
    for (int i = 0; i < N; i++) {
        free(A[i]);
        free(B[i]);
        free(C[i]);
    }

    free(C_check);
    free(A);
    free(B);
    free(C);

    printf("elasped time: %f\n", end-start);

    return 0;
}