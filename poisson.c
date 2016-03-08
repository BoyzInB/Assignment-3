/*Poisson solver using finite differences SOR
 with openmp
 */

/*hai, so much code, such pretty */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <sys/time.h>

#define MAX_SWEEPS	1000000

#define PI	3.1415926536
#define DOUBLE		((unsigned) sizeof(double))
#define ERROR	(-1)
#define sq(x)	( (x)*(x) )
#define EPSILON_M	2.2204460492503131e-16

typedef struct {
    int N;
    int Nx,Ny;
    double dx, i_max;
    double dy, j_max;
    double h;
    double omega, EPSILON;
} GRID;

typedef void	*POINTER;

void implement_BCs(GRID *grid, double *u[]);
void print_solution(GRID *grid, double *u[]);
void SOR(GRID *grid, double *u[], double *f[]);
void fillF(GRID *grid, double *f[]);





int main(void)
{
    struct timeval start, end;
    double diff;
    
    GRID grid;
    double h, mu;
    /* for convergence of iterative methods, EPSILON = epsilon*h^2 */
    double epsilon = 1.e-5;
    int N,i;
    int Nx, Ny;
    
    N = 10;
    grid.Nx = 10;
    grid.Ny = 10;
    grid.N = N;
    int nTemp = N - (N%2 == 0);
    int nxTemp = Nx - (Nx%2 == 0);
    int nyTemp = Ny - (Ny%2 == 0);
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
    printf("SOR omega_opt = %g\n", grid.omega);
    
    printf("hej 1\n");

    /* allocate memory for array u */
    double **u = calloc(Nx,sizeof(*u));
    //double *arr1 = calloc(Ny, sizeof(*arr1));
    double **f = calloc(Nx,sizeof(*f));
    //double *arr2 = calloc(Ny, sizeof(*arr2));
    
    printf("\thej 2\n");

    
  
    
    
    for (i = 0; i < Nx; ++i)
    {
        u[i] = calloc(Ny, sizeof(double));
        f[i] = calloc(Ny, sizeof(double));
    }
    printf("\t\thej 3\n");

    
    //fillF(&grid,f);
    
    //print_solution("f",&grid,f);
    
    
    //fill(&grid, u);
    //    implement_BCs(&grid, u);
    gettimeofday(&start, NULL);
    //SOR(&grid, u,f);
    gettimeofday(&end, NULL);
    
    
    diff = ((end.tv_sec * 1000000 + end.tv_usec)
            - (start.tv_sec * 1000000 + start.tv_usec))/1000000.0;
    printf("Time: %lf sec.\n", diff);
    
    print_solution(&grid, u);
    
    
    
    //free(arr1);
    //free(arr2);
    free(u);
    free(f);
    return 0;
}



void fillF(GRID *grid, double *f[]){
    double dx = grid->dx;
    double dy = grid->dy;
    int N = grid->N;
    int Nx = grid->Nx;
    int Ny = grid->Ny;

    int i,j;
    
    for (i = 0; i<Nx; i++)
        for(j = 0;j<Ny;j++){
            f[i][j] = sin(2*PI*(i-0.5)*dx);
            // printf("x: %f\n",(i-0.5)*dx);
        }
    
}




void implement_BCs(GRID *grid, double *u[])
{
    int N = grid->N, i, j;
    int Nx = grid->Nx;
    int Ny = grid->Ny;


    printf("hej 1\n");

    for (i = 1; i < Nx-1 ; i++) {
        u[0][i] = u[1][i];

        u[Ny-1][i] = u[Ny-2][i];
    }
    printf("\thej 2\n");

    for (j = 1; j < Ny-1; j++) {
        u[j][0] = u[j][1];
        u[j][Nx-1] = u[j][Nx-2];
    }
    u[0][0] = u[1][1];u[Ny-1][0]=u[Ny-2][1];
    u[0][Nx-1] = u[1][Nx-2];u[Ny-1][Nx-1]=u[Ny-2][Nx-2];

     }


void print_solution(GRID *grid, double *u[])
{
    int N = grid->N, i, j;
    int Nx = grid->Nx;
    int Ny = grid->Ny;

    double x, y;
    
    // x, y, u
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 2 ; j++) {
            printf("%.5f ", u[j][i]);
        }
        printf("\n");
    }
    
}

void SOR(GRID *grid, double *u[], double *f[])
{
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
    
    
    printf("hej 1 \n");

    

    double **residual = calloc(Ny,sizeof(*residual));
    double *arr3 = calloc(Nx, sizeof*arr3);
    
    printf("\thej 2\n");

    

    for (i = 0; i < Ny; ++i)
    {
        residual[i] = &arr3[i * (Nx)];
    }
    
    printf("\t\thej 3 \n");

  
    int count;
    
    //Step 1: enforce BC
    printf("\t\t\thej 4 \n");

    implement_BCs(grid,u);
    printf("\t\t\t\thej 5 \n");

#pragma omp parallel
    {
        count = 0;
        while(count < 50000 && norm_residual > EPSILON){

#pragma omp single
            count++;
            //Step 2: compute new u
            
            //Red

#pragma omp for private(j)
            for (i = 1; i < Nx-1; i++){
                for (j = 1 + (i%2!=0); j < Ny-1; j+=2){

                    u[j][i] =(1-omega)*u[j][i];

                    u[j][i] += omega/(2/(dx*dx)+2/(dy*dy));

                    u[i][j] *= ( (u[j+1][i]+u[j-1][i])/(dx*dx)+
                         (u[j][i+1]+u[j][i-1])/(dy*dy) -f[j][i]);
                }
            }
            
            //Black
#pragma omp for private(j)
            for (i = 1; i < Nx-1; i++)
                    for (j = 1 + (i%2==0); j < Ny-1; j+=2)
                        u[j][i] =(1-omega)*u[j][i] + omega/(2/(dx*dx)+2/(dy*dy))*
                        ( (u[j+1][i]+u[j-1][i])/(dx*dx)+
                         (u[j][i+1]+u[j][i-1])/(dy*dy) -f[j][i]);
            
            //Step 3: enforce BC
#pragma omp single nowait
            implement_BCs(grid,u);
            //Step 4: compute the residual
            
            sum = 0.0;
#pragma omp for reduction(+:sum) private(j)
            for (i = 1; i < Nx-1; i++) {
                for (j = 1; j < Ny-1; j++) {
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
    printf("sor done ...\n");
    
    
}
