#include "graph.h"

#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>

#define ISIZE 16
#define STRING_LENGTH 511

/**
 * Graph data structure:
 *
 *   - an array of edges (each edge actually appears twice, once in each
 * direction)
 *
 * @todo fix this so that each edge appears only once
 *
 *   - an array of vertices
 *
 *   - an edge knows its two endpoints and its weight
 *
 *   - a vertex knows a list of its incident edges (each edge is represented
       as an index into the edge array)
 */

GRAPH * new_graph(vertex n, edge m)
{
    GRAPH * g = (GRAPH *)malloc(sizeof(GRAPH));
    g->num_verts = n;
    g->num_edges = 0;
    g->capacity = m;
    g->vertex_array = (VERTEX *)malloc( (n + 1) * sizeof(VERTEX) );
    g->edge_array = (EDGE *) malloc( (m + 1) * sizeof(EDGE) );
    for (vertex i = 1; i <= n; i++) {
      g->vertex_array[i].num_incident_edges = 0;
      g->vertex_array[i].edge_array = (edge *) malloc((ISIZE + 1) * sizeof(edge));
      g->vertex_array[i].x_coord = -1;
      g->vertex_array[i].y_coord = -1;
    }
    /*
      allocate g->edge_array (an array of edge structs)
     */
    return g;
}

// free all data structures
void delete(GRAPH * g )
{
    for (vertex i = 1; i <= g->num_verts; i++)
    {
        free(g->vertex_array[i].edge_array);
    }
    free(g->vertex_array);
    free(g->edge_array);
    free(g);
}

/**
 * Here's where we fix the code to store each edge once.
 */

void add_edge(GRAPH * g, vertex v_one, vertex v_two, edge_weight wt)
{
    vertex u = v_one < v_two ? v_one : v_two, v = v_one > v_two ? v_one : v_two;
    edge new_edge = ++ g->num_edges;
    // size of edge array may have to increase
    if ( new_edge >= g->capacity )
    {
        g->edge_array = realloc(g->edge_array, (g->capacity + ISIZE + 1) * sizeof (EDGE));
        g->capacity += ISIZE;
    }
    g->edge_array[new_edge].vertex_one = u; g->edge_array[new_edge].vertex_two = v; g->edge_array[new_edge].weight = wt;
    vertex j = ++ g->vertex_array[u].num_incident_edges;
    if (j % ISIZE == 0)
        g->vertex_array[u].edge_array = realloc(g->vertex_array[u].edge_array, (j + ISIZE + 1) * sizeof(edge));
    g->vertex_array[u].edge_array[j] = new_edge;
    j = ++ g->vertex_array[v].num_incident_edges;
    if (j % ISIZE == 0)
        g->vertex_array[v].edge_array = realloc(g->vertex_array[v].edge_array, (j + ISIZE + 1) * sizeof(edge));
    g->vertex_array[v].edge_array[j] = new_edge;
}

  // a macro INCIDENT_EDGE( v, i ) to retrieve the i-th incident edge of v

/**
 * reads a graph from a .gph file -- code for this is already in each of the
 * algorithm implementations
 *
 * @todo consolidate everything to use the graph module
 * @todo read and save comments and then write them as well
 */
GRAPH * read_graph( FILE * input_file ) {
  /**
   * @todo fix this so it pays attention to and stores comments before the g
   * part;
   */
  size_t n = 0;
  size_t m = 0;

  char line[STRING_LENGTH];
  /* Scan through comments to get to the "g n m" line; here any string
     that does not begin with 'g' is treated as a comment, i.e., no error
     checking */
  while ( fgets(line, sizeof(line), input_file) ) {
#if defined(DEBUG)
    printf( "reading input, line = %s\n", line );
#endif
    if ( line[0] == 'g' ) {
      sscanf(line, "g %zu %zu", &n, &m);
      break;
    }
  }

  GRAPH * g = new_graph( n, m );

  // now read information about vertices and edges:
  //     e v_one v_two weight     - for an edge
  //     n id x y                    - for vertex id in position (x,y)
  // any other information is simply ignored
  while ( fgets(line, sizeof(line), input_file) ) {
    if ( line[0] == 'e' ) {
      vertex v_one, v_two;
      edge_weight wt;
      sscanf( line, "e "VERTEX_FORMAT" "VERTEX_FORMAT" "EDGE_WEIGHT_FORMAT, &v_one, &v_two, &wt );
      add_edge( g, v_one, v_two , wt );
    }
    else if ( line[0] == 'n' ) {
      vertex id;
      edge_weight x;
      edge_weight y;
      sscanf( line, "n "VERTEX_FORMAT" "EDGE_WEIGHT_FORMAT" "EDGE_WEIGHT_FORMAT,
              &id, &x, &y );
      g->vertex_array[id].x_coord = x;
      g->vertex_array[id].y_coord = y;
    }
  }
  return g;
}

// write a graph to a file in .gph format
void write_graph( FILE * fp, GRAPH * g)
{
  fprintf( fp, "g %zu %zu\n", g->num_verts, g->num_edges);
  for ( vertex v = 1; v <= g->num_verts; v++ ) {
    if ( g->vertex_array[v].x_coord >= 0 ) {
      fprintf( fp, "n "VERTEX_FORMAT" "EDGE_WEIGHT_FORMAT" "EDGE_WEIGHT_FORMAT" \n",
               v, g->vertex_array[v].x_coord, g->vertex_array[v].y_coord );
    }
  }
  for ( edge e = 1; e <= g->num_edges; e++ ) {
    fprintf(fp, "e "VERTEX_FORMAT" "VERTEX_FORMAT" "EDGE_WEIGHT_FORMAT"\n",
            g->edge_array[e].vertex_one,
            g->edge_array[e].vertex_two,
            g->edge_array[e].weight );
  }
}
