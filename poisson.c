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

#include <cdecs.h>

#define MAX_SWEEPS	1000000

#define SOR_METHOD		1
#define GAUSS_SEIDEL_METHOD	2
#define JACOBI_METHOD		3
#define FULL_MATRIX		4
#define BAND_MATRIX		5

typedef struct {
    int N;
    double dx, xmin, xmax;
    double dy, ymin, ymax;
    double h;
    double omega, EPSILON;
} GRID;

void implement_BCs(GRID *grid, double *u[]);
void print_solution(char *string, GRID *grid, double *u[]);
void SOR(GRID *grid, double *u[]);
void jacobi(GRID *grid, double *u[]);
void full_matrix(GRID *grid, double *u[]);
void band_matrix(GRID *grid, double *u[]);
void full_solve(int N, double *A[], double b[]);

int main(void)
{
    double **u;
    GRID grid;
    double xmin = 0., xmax = 1.;
    double ymin = 0., ymax = 1.;
    double h, mu;
    /* for convergence of iterative methods, EPSILON = epsilon*h^2 */
    double epsilon = 1.e-5;
    int N;
    int method = ERROR;
    PFV elliptic;
    char c[100];
    
    printf("Laplace equation solution\n");
    
    fprintf(stderr, "Choose matrix solve method.  Current choices are\n");
    fprintf(stderr, "\tJacobi (J) iteration,\n");
    fprintf(stderr, "\tGauss-Seidel (GS) iteration,\n");
    fprintf(stderr, "\tSOR (SOR) iteration\n");
    fprintf(stderr, "\tfull-matrix (FULL) direct solve\n");
    fprintf(stderr, "\tor band-matrix (BAND) direct solve.\n");
    fprintf(stderr, "Enter choice: ");
    scanf("%s", c);
    if (strcmp(c, "J") == 0 || strcmp(c, "j") == 0) {
        method = JACOBI_METHOD;
        printf("matrix solve method = Jacobi\n");
    }
    else if (strcmp(c, "GS") == 0 || strcmp(c, "gs") == 0) {
        method = GAUSS_SEIDEL_METHOD;
        printf("matrix solve method = Gauss-Seidel\n");
    }
    else if (strcmp(c, "SOR") == 0 || strcmp(c, "sor") == 0) {
        method = SOR_METHOD;
        printf("matrix solve method = SOR\n");
    }
    else if (strcmp(c, "FULL") == 0 || strcmp(c, "full") == 0) {
        method = FULL_MATRIX;
        printf("direct solve method = FULL matrix\n");
    }
    else if (strcmp(c, "BAND") == 0 || strcmp(c, "band") == 0) {
        method = BAND_MATRIX;
        printf("direct solve method = BAND matrix\n");
    }
    else {
        fprintf(stderr, "ERROR: unrecognized method %s\n", c);
        exit(ERROR);
    }
    
    fprintf(stderr, "Enter the number of dx (= number of dy): ");
    scanf("%d", &N);
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
    
    switch (method) {
        case SOR_METHOD:
            elliptic = SOR;
            /* calculate omega_opt for SOR */
            mu = cos(PI*h); /* Jacobi spectral radius */
            grid.omega = 2.*(1.-sqrt(1.-sq(mu)))/sq(mu);
            printf("EPSILON = epsilon*h^2 = %g\n", grid.EPSILON);
            printf("\twhere epsilon = %g\n", epsilon);
            printf("SOR omega_opt = %g\n", grid.omega);
            break;
        case GAUSS_SEIDEL_METHOD:
            elliptic = SOR;
            grid.omega = 1.;
            printf("EPSILON = epsilon*h^2 = %g\n", grid.EPSILON);
            printf("\twhere epsilon = %g\n", epsilon);
            break;
        case JACOBI_METHOD:
            elliptic = jacobi;
            printf("EPSILON = epsilon*h^2 = %g\n", grid.EPSILON);
            printf("\twhere epsilon = %g\n", epsilon);
            break;
        case FULL_MATRIX:
            elliptic = full_matrix;
            break;
        case BAND_MATRIX:
            elliptic = band_matrix;
            break;
        default:
            printf("ERROR: unrecognized method\n");
            exit(ERROR);
            break;
    }
    
    /* allocate memory for array u */
    u = (double **) alloc_matrix(N+1, N+1, DOUBLE);
    zero_matrix(u, N+1, N+1);
    
    implement_BCs(&grid, u);
    (*elliptic)(&grid, u);
    print_solution("solution", &grid, u);
    cpu_time();
    
    return 0;
}

void implement_BCs(GRID *grid, double *u[])
{
    int N = grid->N, i, j;
    
    for (j = 0; j <= N; j++) {
        u[j][0] = u[j][1];
        u[j][N] = u[j][N-1];
    }
    for (i = 1; i < N; i++) {
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
    FILE *file1;
    
    file1 = fopen("laplace.sol", "w");
    
    // x, y, u
    for (i = 0; i <= N; i++) {
        for (j = 0; j <= N; j++) {
            x = grid->xmin + i*grid->dx;
            y = grid->ymin + j*grid->dy;
            fprintf(file1, "%-12g %-12g %-12g\n", x, y, u[j][i]);
        }
    }
    
    fclose(file1);
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

void jacobi(GRID *grid, double *u[])
{
    static double **new_u;
    int N = grid->N, i, j, sweep;
    double EPSILON = grid->EPSILON;
    double sum, residual, norm_residual0, norm_residual = 1.;
    
    new_u = (double **) alloc_matrix(N+1, N+1, DOUBLE);
    for (i = 0; i <= N; i++)
        for (j = 0; j <= N; j++) new_u[j][i] = u[j][i];
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
                new_u[j][i] += residual/4.;
            }
        }
        norm_residual = sum/sq(N-1);
        for (i = 1; i < N; i++)
            for (j = 1; j < N; j++) u[j][i] = new_u[j][i];
    }
    
    printf("number of sweeps = %d\n", sweep);
    printf("norm_residual/norm_residual0 = %g\n",
           norm_residual/norm_residual0);
    if (norm_residual > EPSILON*norm_residual0)
        printf("WARNING: iterative method failed to converge\n");
}

void full_matrix(GRID *grid, double *u[])
{
    int N = grid->N, i, j;
    static double **A, *b;
    int n, modes;
    
    /* allocate memory for arrays */
    modes = sq(N-1); /* unknowns = interior u_ij values */
    A = (double **) alloc_matrix(modes, modes, DOUBLE);
    b = (double *) alloc_vector(modes, DOUBLE);
    
    zero_matrix(A, modes, modes);
    zero_vector(b, modes);
    
    /*
     *	set up A & b in Laplace equation
     *		A u = b
     */
    for (i = 0; i <= N-2; i++) {
        for (j = 0; j <= N-2; j++) {
            n = j*(N-1) + i;
            A[n][n] = -4.;
            if (i > 0) A[n][n-1] = 1.;
            if (i < N-2) A[n][n+1] = 1.;
            if (n >= N-1) A[n][n-(N-1)] = 1.;
            if (n <= sq(N-1)-N) A[n][n+(N-1)] = 1.;
        }
    }
    for (j = 0; j <= N-2; j++) {
        n = j*(N-1);
        b[n] -= u[j+1][0]; /* West BC */
        n = j*(N-1) + (N-2);
        b[n] -= u[j+1][N]; /* East BC */
    }
    for (i = 0; i <= N-2; i++) {
        n = i;
        b[n] -= u[0][i+1]; /* South BC */
        n = (N-2)*(N-1) + i;
        b[n] -= u[N][i+1]; /* North BC */
    }
    
    full_solve(modes, A, b);
    for (i = 0; i <= N-2; i++) {
        for (j = 0; j <= N-2; j++) {
            n = j*(N-1) + i;
            u[j+1][i+1] = b[n];
        }
    }
}

void band_matrix(GRID *grid, double *u[])
{
    printf("band matrix solve not implemented in this version!\n");
    exit(ERROR);
}

void full_solve(int N, double *A[], double b[])
{
    int i, j, k, pivot_row;
    double pivot, multiplier, *rowswap, swap, sum;
    
    /* transform A to upper triangular form */
    /* loop over columns of A */
    for (k = 0; k < N-1; k++) {
        pivot = A[k][k];
        pivot_row = k;
        /* loop over rows of A to find pivot */
        for (i = k+1; i < N; i++) {
            if (ABS(pivot) < ABS(A[i][k])) {
                pivot = A[i][k];
                pivot_row = i;
            }
        }
        if (ABS(pivot) < EPSILON_M) {
            printf("ERROR in solve(): pivot too small\n");
            exit(ERROR);
        }
        if (pivot_row != k) {
            /* swap rows k & pivot_row of A */
            rowswap = A[k];
            A[k] = A[pivot_row];
            A[pivot_row] = rowswap;
            /* swap rows of b */
            swap = b[k];
            b[k] = b[pivot_row];
            b[pivot_row] = swap;
        }
        /* loop over rows of A */
        for (i = k+1; i < N; i++) {
            /* compute multiplier */
            multiplier = A[i][k] = A[i][k]/A[k][k];
            /* compute new row elements of A */
            for (j = k+1; j < N; j++)
                A[i][j] -= multiplier*A[k][j];
            /* compute new elements of b */
            b[i] -= multiplier*b[k];
        }
    }
    
    /* back substitution */
    /* loop over rows in reverse order */
    for (i = N-1; i >= 0; i--) {
        sum = 0.;
        for (j = i+1; j < N; j++)
            sum += A[i][j]*b[j];
        sum = b[i] - sum;
        if (ABS(A[i][i]) > EPSILON_M) b[i] = sum/A[i][i];
        else if (ABS(sum) < EPSILON_M) {
            b[i] = 0.;
            printf("WARNING: sum may be too small in solve()\n");
        }
        else {
            printf("ERROR in solve()\n");
            exit(ERROR);
        }
    }
}