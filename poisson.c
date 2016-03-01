/*
 *				solvers: main.c
 *
 *	This program solves Laplace's equation on the unit square in
 *	finite-difference form using Jacobi, Gauss-Seidel, or SOR iteration,
 *	or a full-matrix or banded-matrix direct solve.
 *
 *	The BCs are Dirichlet with u(B) = 0 except u(x = 1, y) = 1.
 *
 *	In u(x_i,y_j) = u_ij = u[j][i], x[i] runs from xmin to xmax & y[j] runs
 *	from ymin to ymax:
 *		x = xmin + i*dx;
 *		y = ymin + j*dy;
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SWEEPS	1000000

#define PI	3.142
#define DOUBLE		((unsigned) sizeof(double))
#define ERROR	(-1)
#define sq(x)	( (x)*(x) )
#define EPSILON_M	2.2204460492503131e-16

typedef struct {
    int N;
    double dx, xmin, xmax;
    double dy, ymin, ymax;
    double h;
    double omega, EPSILON;
} GRID;

typedef void	*POINTER;

void implement_BCs(GRID *grid, double *u[]);
void print_solution(char *string, GRID *grid, double *u[]);
void SOR(GRID *grid, double *u[]);

POINTER alloc_matrix(int N_rows, int N_columns, unsigned element_size);
POINTER alloc_vector(int N, unsigned element_size);
POINTER Alloc(unsigned N_bytes);
void zero_matrix(double *A[], int N_rows, int N_columns);
void zero_vector(double v[], int N);

int main(void)
{
    GRID grid;
    double xmin = 0., xmax = 1.;
    double ymin = 0., ymax = 1.;
    double h, mu;
    /* for convergence of iterative methods, EPSILON = epsilon*h^2 */
    double epsilon = 1.e-5;
    int N,i;
    void *elliptic;
    char c[100];
    
    N = 10;
    printf("number of dx = number of dy = %d\n", N);
    grid.N = N;
    grid.dx = (xmax-xmin)/N;
    grid.dy = (ymax-ymin)/N;
    h = grid.h = grid.dx;
    grid.xmin = xmin;
    grid.xmax = xmax;
    grid.ymin = ymin;
    grid.ymax = ymax;
    grid.EPSILON = epsilon*sq(h);
    
    elliptic = SOR;
    /* calculate omega_opt for SOR */
    //mu = cos(PI*h); /* Jacobi spectral radius */
    grid.omega = 1.7;
    printf("EPSILON = epsilon*h^2 = %g\n", grid.EPSILON);
    printf("\twhere epsilon = %g\n", epsilon);
    printf("SOR omega_opt = %g\n", grid.omega);
    
    /* allocate memory for array u */
    double **u = calloc((N+1),sizeof(*u));
    double *arr = calloc((N+1)*(N+1), sizeof*arr);
    
    for (i = 0; i < (N+1); ++i)
    {
        u[i] = &arr[i * (N+1)];
    }
    
    implement_BCs(&grid, u);
    SOR(&grid, u);
    print_solution("solution", &grid, u);
    
    free(arr);
    free(u);
    return 0;
}

void implement_BCs(GRID *grid, double *u[])
{
    int N = grid->N, i, j;
    
    for (j = 0; j <= N; j++) {
        u[j][0] = 0.;
        u[j][N] = 1.;
    }
    for (i = 1; i < N; i++) {
        u[0][i] = 0.;
        u[N][i] = 0.;
    }
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


void print_solution(char *string, GRID *grid, double *u[])
{
    int N = grid->N, i, j;
    double x, y;
    
    // x, y, u
    for (i = 0; i <= N; i++) {
        for (j = 0; j <= N; j++) {
            printf("%.2f ", u[j][i]);
        }
        printf("\n");
    }
    
}

void SOR(GRID *grid, double *u[])
{
    Grid newGrid = grid;
    int N = grid->N, i, j, sweep;
    double dx = grid->dx;
    double dy = grid->dy;
    double EPSILON = grid->EPSILON;
    double omega = grid->omega;
    double sum, temp, norm_residual0, norm_residual = 1.;
    
    
    double **residual = calloc((N+1),sizeof(*residual));
    double *arr = calloc((N+1)*(N+1), sizeof*arr);
    
    for (i = 0; i < (N+1); ++i)
    {
        residual[i] = &arr[i * (N+1)];
    }
    
    while(norm_residual > EPSILON){
        //Step 1: enforce BC
        
        //Step 2: compute new u
        for (i = 1; i < N; i++) {
            for (j = 1; j < N; j++) {
                u[j][i] =(1-omega)*u[j][i] + omega/(2/(dx*dx)+2/(dy*dy))*
                ( (u[j+1][i]+u[j-1][i])/(dx*dx)+
                 (u[j][i+1]+u[j][i-1])/(dy*dy) );
            }
        }
        //Step 3: enforce BC
        
        //Step 4: compute the residual
        for (i = 1; i < N; i++) {
            for (j = 1; j < N; j++) {
                residual[j][i] = -((u[j+1][i] - 2*u[j][i] + u[j-1][i])/(dx*dx)
                                   + (u[j][i+1] - 2*u[j][i] + u[j][i-1])/(dy*dy));
            }
        }

        print_solution("residual", &newGrid, residual);

        
        //Step 5: compute it's L2norm
        sum = 0;
        for (i = 1; i < N; i++) {
            for (j = 1; j < N; j++) {
                sum += residual[j][i]*residual[j][i];
            }
        }
        norm_residual = sqrt(1/(N*N)*sum);
        
        printf("%f\n",norm_residual);
    }
    printf("sor done ...\n");
    
    
}
