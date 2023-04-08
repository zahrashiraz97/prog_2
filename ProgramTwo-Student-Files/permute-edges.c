/**
 * @file permute-edges.c
 * @brief a simple filter that reads a graph from standard input, permutes
 * the edges using a given random seed and writes the result to standard
 * output
 * @author Matt Stallmann, 2010/08/09
 * $Id: permute-edges.c 300 2013-11-15 15:36:32Z mfms $
 *
 * @todo Add information in comments that are parsed by read_graph, output
 * by write_graph, and accessed by appropriate functions - see src-new for
 * some indication.
 *   - c # type n m seed
 *   - c | permutation_seed
 *   - c = max_weight
 */

#include "graph.h"
#include "random.h"
#include "common.h"

#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>

/**
 * PRE: A is an array of 'length' items, each of which is 'element_size'
 *      bytes long.
 * POST: A has been randomly permuted
 */
static void permute_edges( GRAPH * g ) {
  EDGE * temp = (EDGE *) malloc( sizeof(EDGE) );
  /** @todo Edge array should probably be indexed starting at 0  */
  for ( edge i = g->num_edges; i > 1; i-- ) {
    /* swap A[i] with a random element among A[1],...,A[i] */
    edge j = genrand_int32() % i + 1;
    if ( j < i ) {
      * temp = g->edge_array[i];
      g->edge_array[i] = g->edge_array[j];
      g->edge_array[j] = * temp;
    }
  }
  free( temp );
}

int main( int argc, char * argv[] ) {
  if ( argc != 2 ) {
    printf( "Usage: %s seed < input > output\n", argv[0] );
    exit( 1 );
  }

  unsigned long seed = strtoul( argv[1], NULL, 10 );
  init_genrand( seed );  

  // save first comment (info about graph used by heuristics)
  char first_line[ BUFFER_SIZE ];
  fgets( first_line, BUFFER_SIZE, stdin );

  GRAPH * g = read_graph( stdin );
  permute_edges( g );

  fputs( first_line, stdout );
  write_graph( stdout, g );
  delete( g );
}

/*  [Last modified: 2021 04 27 at 19:34:49 GMT] */
