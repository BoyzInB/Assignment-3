/**********************************************************************
 * assignment 3
 *
 *
 **********************************************************************/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265


void printArray(const double *A,int n){
    //Print array
    for (int i = 0; i<n; i++){
        for(int j=0;j<n;j++)
            printf("%f  ",A[i*n+j]);
        printf("\n");
    }
}

double sinfunc (double x, double y){
    double res = sin(2*PI*x);
    return res;
}

int main(int argc, char *argv[]) {
    
    int i,j,n = 10;
    
    double *A,*f;
    
    //initialize the array
    A = (double *) malloc(n*n*sizeof(double *));
    f = (double *) malloc(n*n*sizeof(double *));

    
    //initial values ...
    for (i = 0; i<n; i++){
        for(j=0;j<n;j++){
            // fill boundaries
            if(i == 0 || i == n-1)
                A[i*n+j] = 0.0;
            else if(j == 0 || j == n-1)
                A[i*n+j] = 0.0;
            // fill inner
            else
                A[i*n+j] = 1.0;
        }
    }
    
    //f-function
    for (i = 0; i<n; i++){
        for(j=0;j<n;j++){
            f[i*n+j] = sinfunc(i,j);
        }
    }
    //Extra value
    A[(3-1)*n+(4-1)] = 5.0;
    
   
    //Print array
    printArray(A,n);
    printf("\n");
    printArray(f,n);
  
    free(A);
    free(f);

return 0;
}

