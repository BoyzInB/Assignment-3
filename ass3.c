//
//  ass3.c
//  
//
//  Created by Viktor Mattsson on 2016-03-01.
//
//
/*
#include <stdio.h>

#define n=10


typedef struct {
    int N;
    double dx, xmin, xmax;
    double dy, ymin, ymax;
    double h;
    double omega, EPSILON;
} GRID;

void SOR(GRID *grid, double *u[])

void printMatrix((int*) A){
    
}

int main{
    /*int** x;
    int* temp;
    int dimension1_max = 1000, dimension2_max = 1000;
    
    x = malloc(dimension1_max * sizeof(int*));
    temp = malloc(dimension1_max * dimension2_max * sizeof(int));
    for (int i = 0; i < dimension1_max; i++) {
        x[i] = temp + (i * dimension2_max);
    }
    int *A;
    A = (int *) malloc(n*n*sizeof(int *));
    
    
    
    
    

    
    /*free(temp);
    free(x);
    free(A);
    return 0;
}

void SOR(GRID *grid, double *u[])
{
    int N = grid->N, i, j, sweep;
    double EPSILON = grid->EPSILON;
    double omega = grid->omega;
    double sum, residual, norm_residual0, norm_residual = 1.;
    
    sum = 0;
    for (i = 1; i < N; i++) {
        for (j = 1; j < N; j++) {
            residual = (-4.*u[j][i]+u[j+1][i]+u[j-1][i]+
                        u[j][i+1]+u[j][i-1]);
            sum += fabs(residual);
        }
    }
    norm_residual0 = sum/sq(N-1);
    
    for (sweep = 0; sweep<MAX_SWEEPS &&
         norm_residual>EPSILON*norm_residual0; sweep++) {
        sum = 0;
        for (i = 1; i < N; i++) {
            for (j = 1; j < N; j++) {
                residual = (-4.*u[j][i]+u[j+1][i]+u[j-1][i]+
                            u[j][i+1]+u[j][i-1]);
                sum += fabs(residual);
                u[j][i] += omega*residual/4.;
            }
        }
        norm_residual = sum/sq(N-1);
    }
    printf("number of sweeps = %d\n", sweep);
    printf("norm_residual/norm_residual0 = %g\n",
           norm_residual/norm_residual0);
    if (norm_residual > EPSILON*norm_residual0)
     printf("WARNING: iterative method failed to converge\n");
}

*/


/**********************************************************************
 * assignment 3
 *
 *
 **********************************************************************/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define PI 3.14159265


int main(int argc, char *argv[]) {
    
    int i,j,n = 10;
    srand(time(NULL));
    
    double *A;
    A = (double *) malloc(n*n*sizeof(double *));
    
    
    for (i = 0; i<n; i++){
        for(j=0;j<n;j++){
            if(i == 0 || i == n-1)
                A[i*n+j] = 0;
            else if(j == 0 || j == n-1)
                A[i*n+j] = 0;
            else
                A[i*n+j] = 1.0*rand()/RAND_MAX;;
        }
    }
    A[(3-1)*n+(4-1)] = 5;
    
    
    
    printArray(A,n);
    
    impBC(A,n);
    
    printf("\n");
    printArray(A,n);
    
    
    free(A);
    return 0;
}


void printArray(const double *A,int n){
    for (int i = 0; i<n; i++){
        for(int j=0;j<n;j++)
            printf("%.3f  ",A[i*n+j]);
        printf("\n");
    }
}

double sinfunc (double x, double y){
    double res = sin(2*PI*x);
    return res;
}

void impBC (double *A,int n){//Function that implements the boundary conditions
    int i;
    for (i=1; i<n ; i++) {
        A[i*n] = A[i*n+1]; //Column 0
        A[i*n+(n-1)] = A[i*n+(n-2)]; //Column n
        A[i] = A[n+i]; //Row 0
        A[(n-1)*n+i] = A[(n-2)*n+i]; //Row n
    }
    A[0] = A[1]; A[n-1]=A[n-2]; //Corner Points
    A[(n-1)*n] = A[(n-2)*n]; A[n*n-1]= A[n*n-2];
}

void implement_BCs(GRID *grid, double *u[])
{
    int N = grid->N, i, j;
    
    for (j = 0; j <= N; j++) {
        u[j][0] = u[j][1];
        u[j][N] = u[j][N-1];
    }
    for (i = 1; i <= N; i++) {
        u[0][i] = u[1][i];
        u[N][i] = u[N-1][i];
    }
    u[0][0] = u[1][0];u[N][0]=u[N-1][0];
    u[0][N] = u[1][N];u[N][N]=u[N][N-1];
}

void SOR(double *A){
    int i,j;
    for (i=1; i<n; i++) {
        for (j=1; j<n; j++) {
            
        }
    }
}