//
//  ass3.c
//  
//
//  Created by Viktor Mattsson on 2016-03-01.
//
//

#include <stdio.h>

typedef struct {
    int N;
    double dx, xmin, xmax;
    double dy, ymin, ymax;
    double h;
    double omega, EPSILON;
} GRID;

void SOR(GRID *grid, double *u[])

int main{
    int** x;
    int* temp;
    
    x = malloc(dimension1_max * sizeof(int*));
    temp = malloc(dimension1_max * dimension2_max * sizeof(int));
    for (int i = 0; i < dimension1_max; i++) {
        x[i] = temp + (i * dimension2_max);
    }
    
    
    
    
    
    
    
    free(temp);
    free(x);
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
