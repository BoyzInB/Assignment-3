/**********************************************************************
 * assignment 3
 *
 *
 **********************************************************************/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
int n;

void printArray(int A[][]){
    for (int i = 0; i<n; i++){
        for(int j=0;j<n;j++)
            printf("%d",A[i][j]);
        printf("\n");
    }
}
   



int main(int argc, char *argv[]) {
    
    int i,j,n = 10;
    int *A;
    
    //A = (int **)malloc(n*sizeof(int *));
    A = (int *) malloc(n*n*sizeof(int *));
    
    
    for (i = 0; i<n; i++){
        for(j=0;j<n;j++){
            if(i == 0 || i == n-1)
                A[i][j] = 0;
            else if(j == 0 || j == n-1)
                A[i][j] = 0;
            else
                A[i][j] = 1;
        }
    }
    A[3][4] = 5;
    
    printArray(A);
  



return 0;
}

