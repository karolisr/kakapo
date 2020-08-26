CC=gcc

kakapolib.so: fqlib.o fqgzlib.o
	$(CC) -Wall -Werror -std=c99 -lz -O3 -shared -o kakapolib.so fqlib.o fqgzlib.o

fqlib.o: fqlib.c
	$(CC) -Wall -Werror -std=c99 -lz -O3 -fPIC -c -o fqlib.o fqlib.c

fqgzlib.o: fqgzlib.c
	$(CC) -Wall -Werror -std=c99 -lz -O3 -fPIC -c -o fqgzlib.o fqgzlib.c

fq: fq.c fqlib.o
	$(CC) -Wall -Werror -std=c99 -lz -O3 -o fq fq.c fqlib.o

fqgz: fqgz.c fqgzlib.o
	$(CC) -Wall -Werror -std=c99 -lz -O3 -o fqgz fqgz.c fqgzlib.o

all: kakapolib.so fq fqgz

clean:
	rm -f *.o

deep-clean:
	rm -f *.o *.so fq fqgz
