ARGS=250
CC=gcc
CFLAGS=-O2

#all: v1 v2 v3 v4

#v1: v1.c
#	$(CC) $(CFLAGS) -Wall -o v1 v1.c -lm
	
v2: v2.c
	$(CC) $(CFLAGS) -Wall -o v2 v2.c -lm
	
#v3: v3.c
#	$(CC) $(CFLAGS) -mavx2 -lm -Wall -o v3 v3.c
	
#v4: v4.c
#	$(CC) $(CFLAGS) -fopenmp -lm -Wall -o v4 v4.c

run: v2 #v2 v3 v4
#	./v1 $(ARGS)
	./v2 $(ARGS)
#	./v3 $(ARGS)
#	./v4 $(ARGS)

clean:
	rm -f v2 
#v2 v3 v4