/**
 * @file graph.h
 *
 * @brief Defines a generic edge-weighted graph data structure with
 * sufficient information to be used by all the MST algorithms.
 *
 * @author Matt Stallmann
 * @date 2010/07/09
 *
 * Note: Currently serves as a skeleton for Yu Tang's implementation of this
 * data structure and others.
 *
 * $Id: graph.h 300 2013-11-15 15:36:32Z mfms $
 */

#ifndef _GRAPH_H
#define _GRAPH_H

#include <stdio.h>
#include <stdlib.h>

typedef size_t vertex;
typedef size_t edge;
#define VERTEX_FORMAT "%zu"
#define EDGE_FORMAT "%zu"

typedef long edge_weight;
#define MAX_EDGE_WEIGHT LONG_MAX
#define EDGE_WEIGHT_FORMAT "%ld"

/**
 * an edge has two vertices and a weight; by convention, vertex_one will
 * always have a lower number than vertex_two; this may help some algorithms
 * avoid traversing edges twice
 */
typedef struct edge_struct {
  vertex vertex_one;
  vertex vertex_two;
  edge_weight weight;
} EDGE;

/**
 * a vertex has a list of incident edges; both the list (an array) and the
 * number of edges are stored
 */
typedef struct vertex_struct {
  edge num_incident_edges;
  edge * edge_array;
  // the following might be useful for geometric graphs; they are set to -1
  // initially to indicate their absence
  edge_weight x_coord;
  edge_weight y_coord;
} VERTEX;

/**
 * a graph is stored as a list (array) of vertices; the edges can be
 * retrieved via the incidence lists for the vertices; vertices are numbered
 * from 1 to num_verts - vertex 0 is not a legal vertex (many graph input
 * files give a special meaning to a 0.
 */
typedef struct graph_struct {
  size_t num_verts;
  size_t num_edges;
  size_t capacity;
  VERTEX * vertex_array;
  EDGE * edge_array;
} GRAPH;

/* functions to implement (in graph.c, except for the macro) */

/**
 * allocates all data structures
 */
GRAPH * new_graph( vertex n, edge m ); // n,m = # vertices, edges

/**
 * frees all data structures
 */
void delete( GRAPH * g );

/**
 * adds an edge
 */
void add_edge(GRAPH * g, vertex v_one, vertex v_two, edge_weight wt);

/* #define INCIDENT_EDGE( v, i ) // to retrieve the i-th incident edge of v */

/**
 * @return a pointer to a graph corresponding to the gph formatted
 * file in file stream fp
 */
GRAPH * read_graph( FILE * fp );

/**
 * writes graph g to file steam fp in .gph format
 */
void write_graph( FILE * fp, GRAPH * g );

/**
 * sorts edges of g lexicographically
 */
void sort_edges( GRAPH * g );

/**
 * retrieve the i-th edge incident on vertex v in graph g
 */

#define INCIDENT_EDGE( g, v, i ) g->vertex_array[v].edge_array[i]

#define EDGE_AT( g, i ) g->edge_array[i]

#define OTHER( g, i, u ) ( (u) ^ (g->edge_array[i].vertex_one) ^ (g->edge_array[i].vertex_two) )

#endif

/*  [Last modified: 2021 04 29 at 18:48:29 GMT] */
