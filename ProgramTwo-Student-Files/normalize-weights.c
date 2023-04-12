/**
 * @file normalize-weights.c
 * @brief a simple filter that reads a graph from standard input, gives new
 * weights to the edges so that the range of edge weights is [1,max] and
 * writes the result to standard output. The value of max is a command-line
 * argument.
 *
 * Vertex positions, if they exist, will be normalized as well, using the
 * same ratio as the edge weights.
 *
 * @author Matt Stallmann, 2010/08/09
 * $Id: normalize-weights.c 300 2013-11-15 15:36:32Z mfms $
 */

#include "common.h"
#include "graph.h"
#include "random.h"

#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

/**
 * @return original_value * normalization_factor; both original value and
 * result are integers and there's a random perturbation in case
 * normalization_factor >> 1
 * return value will always be >= 1
 */
static edge_weight normalize(edge_weight original_value, double normalization_factor) {
    double normalized_value = original_value;
    if ( normalization_factor <= 1 ) {
        normalized_value = original_value * normalization_factor;
    }
    else {
        normalized_value = ((double) original_value + genrand_real2())
            * normalization_factor;
    }
    return (edge_weight) fmax(1, floor(normalized_value));
}

/**
 * @param g a graph
 * @return the maximum edge weight or vertex coordinate in g
 */
static edge_weight get_max_value( GRAPH * g ) {
  edge_weight current_max_value = 0;
  for ( size_t i = 1; i <= g->num_edges; i++ ) {
    if ( g->edge_array[i].weight > current_max_value )
      current_max_value = g->edge_array[i].weight;
  }

  for ( size_t v = 1; v <= g->num_verts; v++ ) {
    if ( g->vertex_array[v].x_coord > current_max_value )
      current_max_value = g->vertex_array[v].x_coord;
    if ( g->vertex_array[v].y_coord > current_max_value )
      current_max_value = g->vertex_array[v].y_coord;
  }

  return current_max_value;
}

static void normalize_vertices( GRAPH * g, double normalization_factor ) {
  for ( size_t v = 1; v <= g->num_verts; v++ ) {
    if ( g->vertex_array[v].x_coord >= 0 ) {
        // note: normalize may return a long, while x and y
        // coordinates are ints
        // not a big deal: coordinates are used only in, for example,
        // animations, where they will be modified using the
        // gph2graphml.py script
        g->vertex_array[v].x_coord
            = normalize( g->vertex_array[v].x_coord, normalization_factor );
        g->vertex_array[v].y_coord
            = normalize( g->vertex_array[v].y_coord, normalization_factor );
    }
  }
}

static void normalize_edges( GRAPH * g, double normalization_factor ) {
  for ( size_t i = 1; i <= g->num_edges; i++ ) {
    g->edge_array[i].weight = normalize( g->edge_array[i].weight,
                                         normalization_factor );
  }
}

int main( int argc, char * argv[] ) {
  if ( argc != 2 ) {
    printf( "Usage: %s desired_max_weight < input > output\n", argv[0] );
    exit( 1 );
  }

  edge_weight desired_max_value = strtoul( argv[1], NULL, 10 );
  /* int seed = atoi( argv[2] ); */
  /* init_genrand( seed ); */

  // save first comment (info about graph used by heuristics)
  char first_line[ BUFFER_SIZE ];
  fgets( first_line, BUFFER_SIZE, stdin );

  GRAPH * g = read_graph( stdin );

  edge_weight max_value = get_max_value( g );

  double normalization_factor = desired_max_value / (double) max_value;
  fprintf( stderr,
           " == Normalizing values: actual = "EDGE_WEIGHT_FORMAT
           ", desired = "EDGE_WEIGHT_FORMAT
           ", factor = %e\n",
          max_value, desired_max_value, normalization_factor );

  normalize_vertices( g, normalization_factor );
  normalize_edges( g, normalization_factor );

  fputs( first_line, stdout );
  write_graph( stdout, g );
  delete( g );
}

/*  [Last modified: 2021 04 29 at 19:20:08 GMT] */
