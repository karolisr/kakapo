CC=gcc

# OPTS=-Wall -Werror -std=c99 -O3
OPTS=-std=c99 -O3

DIR_LIB=../lib
DIR_BIN=../bin

DIR_TEST_DATA=../../../../tests/test_data

OS := $(shell uname -s)
ifeq ($(OS),Linux)
	kkparch := linux
endif
ifeq ($(OS),Darwin)
	kkparch := mac
endif
ARCH := $(shell uname -m)
KAKAPOLIB := kakapolib_$(kkparch)_${ARCH}.so

$(info Compiling $(KAKAPOLIB))

all: $(KAKAPOLIB)

$(KAKAPOLIB): gzlib.o fqlib.o
	$(CC) $(OPTS) -lz -shared -o $(KAKAPOLIB) gzlib.o fqlib.o
	@if [ "$(OS)" = "Linux" ]; then ld -lz -shared -o $(KAKAPOLIB) gzlib.o fqlib.o; fi

fqlib.o: fqlib.c
	$(CC) $(OPTS) -fPIC -c -o fqlib.o fqlib.c

gzlib.o: gzlib.c
	$(CC) $(OPTS) -fPIC -c -o gzlib.o gzlib.c

fq_avg_read_len: fq_avg_read_len.c gzlib.o fqlib.o
	$(CC) $(OPTS) -o fq_avg_read_len fq_avg_read_len.c gzlib.o fqlib.o -lz

fq_avg_read_len_gz: fq_avg_read_len_gz.c gzlib.o fqlib.o
	$(CC) $(OPTS) -o fq_avg_read_len_gz fq_avg_read_len_gz.c gzlib.o fqlib.o -lz

test: all fq_avg_read_len fq_avg_read_len_gz
	time ./fq_avg_read_len $(DIR_TEST_DATA)/FASTQ_300000.fq
	time ./fq_avg_read_len_gz $(DIR_TEST_DATA)/FASTQ_300000.fq.gz
	mkdir -p $(DIR_BIN)
	mv fq_avg_read_len fq_avg_read_len_gz $(DIR_BIN)

clean:
	rm -f *.o

install: all clean
	mkdir -p $(DIR_LIB)
	mkdir -p $(DIR_BIN)
	mv *.so $(DIR_LIB)

deep-clean: clean
	rm -f *.so fq_avg_read_len fq_avg_read_len_gz
	rm -rf $(DIR_LIB)
	rm -rf $(DIR_BIN)
