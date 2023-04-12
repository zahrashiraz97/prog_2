## makefile - for minimum spanning tree algorithms and graph generators.
# $Id: Makefile 264 2010-08-10 01:17:36Z mfms $

COMMON_HEADERS = common.h nextgen.h random.h graph.h makefile
CC = gcc
LINKER = gcc

# choose the first set of options if in (serious) debugging mode,
# the second if gdb is to be used and the third in production mode

FLAGS = -std=c99 -g -Wall -W -pedantic
#FLAGS = -std=c99 -Wall -W -O3 -pedantic

# the following allows you to do "make DFLAGS=..." on the command line; for
# example, replace the ... with -DDEBUG
DFLAGS =

EXECUTABLES = nextgen permute-edges normalize-weights triangulation

all: $(EXECUTABLES)

.c.o :
	$(CC) -c $(FLAGS) $(DFLAGS) $< 

nextgen: random_graph.o geometric.o nextgen.o random.o
	$(CC) random_graph.o geometric.o nextgen.o random.o -lm -o nextgen

permute-edges: permute-edges.o graph.o random.o
	$(CC) permute-edges.o graph.o random.o -o permute-edges

normalize-weights: normalize-weights.o graph.o random.o
	$(CC) normalize-weights.o graph.o random.o -lm -o normalize-weights

triangulation: triangulation.o graph.o
	$(CC) triangulation.o graph.o random.o -o triangulation

graph.o: $(COMMON_HEADERS) graph.c

nextgen.o: $(COMMON_HEADERS) nextgen.c

random.o: $(COMMON_HEADERS) random.c

random_graph.o: $(COMMON_HEADERS) random_graph.c

geometric.o: $(COMMON_HEADERS) geometric.c

permute-edges.o: $(COMMON_HEADERS) permute-edges.c

normalize-weights.o: $(COMMON_HEADERS) normalize-weights.c

triangulation.o: $(COMMON_HEADERS) triangulation.c

clean:
	rm -f *.o $(EXECUTABLES)
	rm -rf *.dSYM
