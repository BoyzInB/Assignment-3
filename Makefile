CC=gcc
LIBS=-lm
FLAGS= -O3 -pthread 

puss: poisson.c
	$(CC) $(FLAGS) -o puss poisson.c

clean:
	rm -f *.o *~ puss join puss \
	puss
