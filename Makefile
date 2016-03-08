CC=gcc
LIBS=-lm
FLAGS= -fopenmp 

poiss: poisson.c
	$(CC) $(FLAGS) -o poiss poisson.c

