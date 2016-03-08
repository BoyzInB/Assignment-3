CC=gcc
LIBS=-lm
FLAGS= 

poiss: poisson.c
	$(CC) $(FLAGS) -o poiss poisson.c

