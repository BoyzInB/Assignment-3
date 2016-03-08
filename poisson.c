/*Poisson solver using finite differences SOR
 with openmp
 
 Assignment 3 Programmering av parallelldatorer, VT-16
 Aleksandar Senek
 Viktor Mattsson
 Robin Eriksson
 */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <sys/time.h>

#define PI	3.1415926536

typedef struct {
    int N;
    int Nx,Ny;
    double dx, i_max;
    double dy, j_max;
    double h;
    double omega, EPSILON;
} GRID;

void implement_BCs(GRID *grid, double *u[]);
void print_solution(GRID *grid, double *u[]);
void SOR(GRID *grid, double *u[], double *f[]);
void fillF(GRID *grid, double *f[]);
void fillOnes(GRID *grid, double *f[]);

int main(void)
{//Setting up
    struct timeval start, end;
    double diff;
    
    GRID grid;
    double h, mu;
    double epsilon = 1.e-5;
    int N,j;
    int Nx = 10, Ny = 10;
    
    N = 10;
    grid.Nx = Nx;
    grid.Ny = Ny;
    grid.N = N;
    double a=1.,b=1.;
    grid.i_max = Nx - 2;
    grid.j_max = Ny - 2;
    grid.dx = a/grid.i_max;
    grid.dy = b/grid.j_max;
    
    h = grid.h = grid.dx;
    
    grid.EPSILON = 1e-13;
    
    
    grid.omega = 1.7;
    printf("EPSILON = epsilon*h^2 = %g\n", grid.EPSILON);
    printf("\twhere epsilon = %g\n", epsilon);
    printf("SOR omega = %g\n", grid.omega);
    
    
    /* allocate memory for array u */
    double **u = calloc(Nx,sizeof(*u));
    double **f = calloc(Nx,sizeof(*f));
    
    
    for (j = 0; j < Nx; ++j)
    {
        u[j] = calloc(Ny, sizeof(double));
        f[j] = calloc(Ny, sizeof(double));
    }
    
    /*
     //Task 1
    fillOnes(&grid,u);
    u[2][3] = 5;
    print_solution(&grid, u);
     */
    
    fillF(&grid,f);
    
    gettimeofday(&start, NULL);// For timing
    
    SOR(&grid, u,f); // Solution.
    
    gettimeofday(&end, NULL);
    
    
    diff = ((end.tv_sec * 1000000 + end.tv_usec)
            - (start.tv_sec * 1000000 + start.tv_usec))/1000000.0;
    printf("Time: %lf sec.\n", diff);
    
    print_solution(&grid, u);
    
    
    
    /* Free memory */
    for (int i = 0; i < Nx; i++){
        free( f[i] );
        free( u[i] );
    }
    free(f);
    free(u);
    return 0;
}





void implement_BCs(GRID *grid, double *u[])
{ //Implements the boundary conditions as described
    int N = grid->N, i, j;
    int Nx = grid->Nx;
    int Ny = grid->Ny;
    
    for (i = 1; i < Ny-1 ; i++) {
        u[0][i] = u[1][i];
        u[Nx-1][i] = u[Nx-2][i];
    }
    
    for (j = 1; j < Nx-1; j++) {
        u[j][0] = u[j][1];
        u[j][Ny-1] = u[j][Ny-2];
    }
    u[0][0] = u[1][1];u[Nx-1][0]=u[Nx-2][1];
    u[0][Ny-1] = u[1][Ny-2];u[Nx-1][Ny-1]=u[Nx-2][Ny-2];
    
}


void print_solution(GRID *grid, double *u[])
{ //Prints the 2D array
    int N = grid->N, i, j;
    int Nx = grid->Nx;
    int Ny = grid->Ny;
    
    double x, y;
    printf("\nNx %d, Ny %d\n",Nx,Ny);
    
    for (i = 0; i < Ny; i++) {
        for (j = 0; j < Nx  ; j++) {
            printf("%.5f ", u[j][i]);
        }
        printf("\n");
    }
    
}

void fillOnes(GRID *grid, double *f[])
{ //Fills the 2D array with ones (task 1)
    double dx = grid->dx;
    double dy = grid->dy;
    int N = grid->N;
    int Nx = grid->Nx;
    int Ny = grid->Ny;
    
    int i,j;
    
    for (i = 1; i<Ny-1; i++)
        for(j = 1;j<Nx-1;j++){
            f[j][i] = 1;
        }
}

void fillF(GRID *grid, double *f[])
{ //Fills the 2D array with values from the sin function
    double dx = grid->dx;
    double dy = grid->dy;
    int N = grid->N;
    int Nx = grid->Nx;
    int Ny = grid->Ny;
    
    int i,j;
    
    for (i = 0; i<Ny; i++)
        for(j = 0;j<Nx;j++){
            f[j][i] = sin(2*PI*(j-0.5)*dx);
        }
}

void SOR(GRID *grid, double *u[], double *f[])
{ //Succesisve Over-Relaxation
    int N = grid->N, i, j;
    int Nx = grid->Nx;
    int Ny = grid->Ny;
    
    
    double dx = grid->dx;
    double dy = grid->dy;
    int i_max = grid->i_max;
    int j_max = grid->j_max;
    double EPSILON = grid->EPSILON;
    double omega = grid->omega;
    double sum, temp, norm_residual0, norm_residual = 1.;
    
    double **residual = calloc(Nx,sizeof(*residual));

    for (i = 0; i < Nx; ++i)
    {
        residual[i] = calloc(Ny, sizeof(double));
        
    }
    
    int count;
    
    //Step 1: enforce BC
    
    implement_BCs(grid,u);
    
#pragma omp parallel
    {
        count = 0;
        while(count < 50000 && norm_residual > EPSILON){
            
#pragma omp single
            count++;
            //Step 2: compute new u
            
            //Red
            
#pragma omp for private(j)
            for (i = 1; i < Ny-1; i++)
                for (j = 1 + (i%2!=0); j < Nx-1; j+=2)
                    u[j][i] =(1-omega)*u[j][i] + omega/(2/(dx*dx)+2/(dy*dy))*
                    ( (u[j+1][i]+u[j-1][i])/(dx*dx)+
                     (u[j][i+1]+u[j][i-1])/(dy*dy) -f[j][i]);
            
            //Black
#pragma omp for private(j)
            for (i = 1; i < Ny-1; i++)
                for (j = 1 + (i%2==0); j < Nx-1; j+=2)
                    u[j][i] =(1-omega)*u[j][i] + omega/(2/(dx*dx)+2/(dy*dy))*
                    ( (u[j+1][i]+u[j-1][i])/(dx*dx)+
                     (u[j][i+1]+u[j][i-1])/(dy*dy) -f[j][i]);
            
            //Step 3: enforce BC
#pragma omp single nowait
            implement_BCs(grid,u);
            //Step 4: compute the residual
            
            sum = 0.0;
#pragma omp for reduction(+:sum) private(j)
            for (i = 1; i < Ny-1; i++) {
                for (j = 1; j < Nx-1; j++) {
                    residual[j][i] = f[j][i]-((u[j+1][i] - 2*u[j][i] + u[j-1][i])/(dx*dx)
                                              + (u[j][i+1] - 2*u[j][i] + u[j][i-1])/(dy*dy));
                    sum = residual[j][i]*residual[j][i];
                }
            }
            //Step 5: compute it's L2norm
#pragma omp single
            norm_residual = sqrt(1.0/((double)i_max*(double)j_max)*sum);
            
        }
    }
    
    printf("count: %d\n", count);
    printf("Sor done!\n");
    
    
}