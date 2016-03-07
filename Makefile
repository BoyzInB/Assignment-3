CC=gcc
LIBS=-lm
FLAGS= -O3 -fopenmp 

poiss: poisson.c
	$(CC) $(FLAGS) -o poiss poisson.c

poiss2: poisson.c
	$(CC) -o poiss2 poisson.c
