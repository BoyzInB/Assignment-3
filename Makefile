CC=gcc
LIBS=-lm
FLAGS= -O3 -fopenmp 

poiss: poisson.c
	$(CC) $(FLAGS) -o poiss poisson.c

