/*
*				cdecs.h
*	Copyright 2008 by Carl L. Gardner.  All rights reserved.
*
*	Contains useful extensions to C.
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef void	*POINTER;	/* pointer to an unknown data type */

#define sq(x)	( (x)*(x) )
#define cube(x)	( (x)*(x)*(x) )
#define	cot(x)	( 1./tan(x) )
#define PI	3.14159265358979323846

#define EPSILON_M	2.2204460492503131e-16

#define YES	(1)
#define NO	(0)
#define TRUE	(1)
#define FALSE	(0)
#define ERROR	(-1)

#define ABS(x)		( ((x) >= 0) ? (x) : -(x) )
#define max(a,b)	( ((a) > (b)) ? (a) : (b) )
#define min(a,b)	( ((a) < (b)) ? (a) : (b) )

#define INT		((unsigned) sizeof(int))
#define FLOAT		((unsigned) sizeof(float))
#define DOUBLE		((unsigned) sizeof(double))
#define CHAR		((unsigned) sizeof(char))

typedef void	(*PFV)();	/* pointer to a function returning void */
typedef int	(*PFI)();	/* pointer to a function returning an int */
typedef float	(*PFF)();	/* pointer to a function returning a float */
typedef double	(*PFD)();	/* pointer to a function returning a double */
typedef char	(*PFC)();	/* pointer to a function returning a char */

/* function prototypes in library lib */
/* alloc.c */
POINTER alloc_matrix(int N_rows, int N_columns, unsigned element_size);
POINTER alloc_vector(int N, unsigned element_size);
POINTER Alloc(unsigned N_bytes);
void zero_matrix(double *A[], int N_rows, int N_columns);
void zero_vector(double v[], int N);
/* time.c */
void cpu_time(void);
