/**********************************************************************
 * assignment 3
 *
 *
 **********************************************************************/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>


void printArray(const int *A,int n){
    for (int i = 0; i<n; i++){
        for(int j=0;j<n;j++)
            printf("%d",A[i*n+j]);
        printf("\n");
    }
}



int main(int argc, char *argv[]) {
    
    int i,j,n = 10;
    
    int *A;
    A = (int *) malloc(n*n*sizeof(int *));
    
    
    for (i = 0; i<n; i++){
        for(j=0;j<n;j++){
            if(i == 0 || i == n-1)
                A[i*n+j] = 0;
            else if(j == 0 || j == n-1)
                A[i*n+j] = 0;
            else
                A[i*n+j] = 1;
        }
    }
    A[(3-1)*n+(4-1)] = 5;
    
   
    
    printArray(A,n);
  



return 0;
}

