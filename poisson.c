/*Poisson solver using finite differences SOR
 with openmp
 */

/*hai, so much code, such pretty */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#define MAX_SWEEPS	1000000

#define PI	3.1415926536
#define DOUBLE		((unsigned) sizeof(double))
#define ERROR	(-1)
#define sq(x)	( (x)*(x) )
#define EPSILON_M	2.2204460492503131e-16

typedef struct {
    int N;
    double dx, i_max;
    double dy, j_max;
    double h;
    double omega, EPSILON;
} GRID;

typedef void	*POINTER;

void implement_BCs(GRID *grid, double *u[]);
void print_solution(char *string, GRID *grid, double *u[]);
void SOR(GRID *grid, double *u[], double *f[]);
void fillF(GRID *grid, double *f[]);



int main(void)
{
    GRID grid;
    double h, mu;
    /* for convergence of iterative methods, EPSILON = epsilon*h^2 */
    double epsilon = 1.e-5;
    int N,i;
    char c[100];
    
    N = 10;
    grid.N = N;
    int nTemp = N - (N%2 == 0);
    double a=1.,b=1.;
    grid.i_max = N - 2;
    grid.j_max = N - 2;
    grid.dx = a/grid.i_max;
    grid.dy = b/grid.j_max;
    
    h = grid.h = grid.dx;
    
    grid.EPSILON = 1e-13;
    
    
    grid.omega = 1.7;
    printf("EPSILON = epsilon*h^2 = %g\n", grid.EPSILON);
    printf("\twhere epsilon = %g\n", epsilon);
    printf("SOR omega_opt = %g\n", grid.omega);
    
    /* allocate memory for array u */
    double **u = calloc((N),sizeof(*u));
    double *arr1 = calloc((N)*(N), sizeof*arr1);
    double **f = calloc((N),sizeof(*f));
    double *arr2 = calloc((N)*(N), sizeof*arr2);
    
    
    for (i = 0; i < (N); ++i)
    {
        u[i] = &arr1[i * (N)];
        f[i] = &arr2[i * (N)];
    }
    
    fillF(&grid,f);
    
    //print_solution("f",&grid,f);
    
    
    //fill(&grid, u);
    //    implement_BCs(&grid, u);
    SOR(&grid, u,f);
    print_solution("solution", &grid, u);
    
    
    
    free(arr1);
    free(arr2);
    free(u);
    free(f);
    return 0;
}



void fillF(GRID *grid, double *f[]){
    double dx = grid->dx;
    double dy = grid->dy;
    int N = grid->N;
    
    for (int i = 0; i<N; i++)
        for(int j = 0;j<N;j++){
            f[i][j] = sin(2*PI*(i-0.5)*dx);
            // printf("x: %f\n",(i-0.5)*dx);
        }
    
}




void implement_BCs(GRID *grid, double *u[])
{
    int N = grid->N, i, j;
    
    for (i = 1; i < N-1 ; i++) {
        u[0][i] = u[1][i];
        u[N-1][i] = u[N-2][i];
    }
    for (j = 1; j < N-1; j++) {
        u[j][0] = u[j][1];
        u[j][N-1] = u[j][N-2];
    }
    u[0][0] = u[1][1];u[N-1][0]=u[N-2][1];
    u[0][N-1] = u[1][N-2];u[N-1][N-1]=u[N-2][N-2];
}


void print_solution(char *string, GRID *grid, double *u[])
{
    int N = grid->N, i, j;
    double x, y;
    
    // x, y, u
    for (i = 0; i < N; i++) {
        for (j = 0; j < N ; j++) {
            printf("%.5f ", u[j][i]);
        }
        printf("\n");
    }
    
}

void SOR(GRID *grid, double *u[], double *f[])
{
    int N = grid->N, i, j, sweep;
    double dx = grid->dx;
    double dy = grid->dy;
    int i_max = grid->i_max;
    int j_max = grid->j_max;
    double EPSILON = grid->EPSILON;
    double omega = grid->omega;
    double sum, temp, norm_residual0, norm_residual = 1.;
    
    double **residual = calloc((N),sizeof(*residual));
    double *arr3 = calloc((N)*(N), sizeof*arr3);
    
    for (i = 0; i < (N); ++i)
    {
        residual[i] = &arr3[i * (N)];
    }
    
    
    int count = 0;
    
    
    while(count < 10000 && norm_residual > EPSILON){
        count++;
        //Step 1: enforce BC
        implement_BCs(grid,u);
        //Step 2: compute new u
        
        //Red
        for (i = 1; i < N-1; i+=2) {
            for (j = 1; j < N-1; j+=2) {
                u[j][i] =(1-omega)*u[j][i] + omega/(2/(dx*dx)+2/(dy*dy))*
                ( (u[j+1][i]+u[j-1][i])/(dx*dx)+
                 (u[j][i+1]+u[j][i-1])/(dy*dy) -f[j][i]);
            }
        }
        
        
        for (i = 2; i < N-1; i+=2) {
            for (j = 2; j < N-1; j+=2) {
                u[j][i] =(1-omega)*u[j][i] + omega/(2/(dx*dx)+2/(dy*dy))*
                ( (u[j+1][i]+u[j-1][i])/(dx*dx)+
                 (u[j][i+1]+u[j][i-1])/(dy*dy) -f[j][i]);
            }
        }
        
        //Black
        for (i = 1; i < N-1; i+=2) {
            for (j = 2; j < N-1; j+=2) {
                u[j][i] =(1-omega)*u[j][i] + omega/(2/(dx*dx)+2/(dy*dy))*
                ( (u[j+1][i]+u[j-1][i])/(dx*dx)+
                 (u[j][i+1]+u[j][i-1])/(dy*dy) -f[j][i]);
            }
        }
        
        
        for (i = 2; i < N-1; i+=2) {
            for (j = 1; j < N-1; j+=2) {
                u[j][i] =(1-omega)*u[j][i] + omega/(2/(dx*dx)+2/(dy*dy))*
                ( (u[j+1][i]+u[j-1][i])/(dx*dx)+
                 (u[j][i+1]+u[j][i-1])/(dy*dy) -f[j][i]);
            }
        }
        
        
        
        
        
        
        
        //Step 3: enforce BC
        implement_BCs(grid,u);
        //Step 4: compute the residual
        {
            for (i = 1; i < N-1; i++) {
                for (j = 1; j < N-1; j++) {
                    residual[j][i] = f[j][i]-((u[j+1][i] - 2*u[j][i] + u[j-1][i])/(dx*dx)
                                              + (u[j][i+1] - 2*u[j][i] + u[j][i-1])/(dy*dy));
                }
            }
            
            
            //Step 5: compute it's L2norm
            sum = 0.0;
            for (i = 1; i < N-1; i++) {
                for (j = 1; j < N-1; j++) {
                    sum += residual[j][i]*residual[j][i];
                }
            }
        }
        norm_residual = sqrt(1.0/((double)i_max*(double)j_max)*sum);
        
        
        //print_solution("solution", grid, u);
        
    }
    printf("count: %d\n", count);
    printf("sor done ...\n");
    
    
}
