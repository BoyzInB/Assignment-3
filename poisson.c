/*
 *				solvers: main.c
 *
 *	This program solves Laplace's equation on the unit square in
 *	finite-difference form SOR iteration,
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
void fill(GRID *grid, double *u[]);
void fillFunc(GRID *grid, double *u[]);
double sinfunc (double x, double y);

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
    fillFunc(&grid,u);
    print_solution("solution", &grid, u);
    implement_BCs(&grid, u);
    SOR(&grid, u);
    print_solution("solution", &grid, u);
    
    free(arr);
    free(u);
    return 0;
}

/*void implement_BCs(GRID *grid, double *u[])
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
 }*/

void fill(GRID *grid, double *u[]){
    int N = grid->N,i,j;
    
    for (i=N/2-1;i<N/2+2;i++)
        for (j=N/2-1;j<N/2+2; j++)
            u[i][j] = 10.;
}

void fillFunc(GRID *grid, double *u[]){
    int N = grid->N,i,j;
    double dx = grid->dx, dy=grid->dy;
    for (i=0;i<N;i++)
        for (j=0;j<N; j++)
            u[i][j] = sinfunc(dx*i,dy*j);
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

double sinfunc (double x, double y){
    double res = sin(2*PI*x);
    return res;
}

void SOR(GRID *grid, double *u[])
{
    int N = grid->N, i, j, sweep;
    double dx = grid->dx;
    double dy = grid->dy;
    double EPSILON = grid->EPSILON;
    double omega = grid->omega;
    double sum, temp, norm_residual0, norm_residual = 1.;
    int count=0;
    
    double **residual = calloc((N+1),sizeof(*residual));
    double *arr = calloc((N+1)*(N+1), sizeof*arr);
    
    for (i = 0; i < (N+1); ++i)
    {
        residual[i] = &arr[i * (N+1)];
    }
    
    while(norm_residual > EPSILON){
        //Step 1: enforce BC
        implement_BCs(grid,u);
        //Step 2: compute new u
        for (i = 1; i < N; i++) {
            for (j = 1; j < N; j++) {
                u[i][j] =(1-omega)*u[i][j] + omega/(2/(dx*dx)+2/(dy*dy))*
                ( (u[i+1][j]+u[i-1][j])/(dx*dx)+
                 (u[i][j+1]+u[i][j-1])/(dy*dy) );
            }
        }
        //Step 3: enforce BC
        implement_BCs(grid,u);
        //Step 4: compute the residual
        for (i = 1; i < N; i++) {
            for (j = 1; j < N; j++) {
                residual[i][j] = -((u[i+1][j] - 2*u[i][j] + u[i-1][j])/(dx*dx)
                                   + (u[i][j+1] - 2*u[i][j] + u[i][j-1])/(dy*dy));
            }
        }
        
        //print_solution("residual", grid, residual);
        
        
        
        //Step 5: compute it's L2norm
        sum = 0;
        for (i = 1; i < N; i++) {
            for (j = 1; j < N; j++) {
                sum += residual[j][i]*residual[j][i];
            }
        }
        printf("%.7f\n",sum);
        norm_residual = sqrt(1.0/((double)N*(double)N)*sum);
        
        printf("%.7f\n",norm_residual);
        count++;
        printf("Count: %d\n",count);
    }
    printf("sor done ...\n");
    
    
}