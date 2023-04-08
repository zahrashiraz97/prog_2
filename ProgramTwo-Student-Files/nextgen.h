/*******************************************************************
 * Barry Peddycord
 *
 * 	Header file for Generator #2
 ******************************************************************/

#ifndef __NEXTGEN_H
#define __NEXTGEN_H

/**
 * @todo the following also appear in graph.h; eventually both
 * random_graph.c and geometric.c should use those definitions
 */

typedef size_t vertex;
typedef size_t edge;
#define VERTEX_FORMAT "%zu"
#define EDGE_FORMAT "%zu"
typedef long edge_weight;
#define MAX_EDGE_WEIGHT LONG_MAX
#define EDGE_WEIGHT_FORMAT "%ld"

/** Name of output file: records information about algorithm run if random
    graph is fed directly to an algorithm - DIRECT_INPUT - or stores the
    graph in a simple format */
extern char * ofname;

/** actual number of edges (in case this is a geometric graph) */
extern size_t actual_number_of_edges;

/** Used only to print a comment at the top of the file for a randomly
    generated graph */
extern uint seed;

#ifndef DIRECT_INPUT
/** File for storing the graph */
extern FILE * outfile;
#endif

/**
 * @brief Generates a random graph with randomly generated and connected
 * vertices
 * @param num_verts number of vertices
 * @param num_edges approximate number of edges (no guarantee) 
 */
void random_graph( size_t num_verts, size_t num_edges );

/**
 * @brief Generates a randomly generated geometric graph; the graph is
 * connect if it's sufficiently dense
 * @param num_verts number of vertices
 * @param num_edges approximate number of edges (no guarantee) 
 * @return true if the graph is connected
 */
bool geometric_graph( size_t num_verts, size_t num_edges );

/**
 * @brief initialize the memory needed for Boruvka's algorithm.
 * @param n Number of vertices.
 * @param m Number of edges.
 */
void init_graph( size_t n, size_t m );

/**
 *	@brief Adds an edge to the graph
 *	@param v_one the vertex at one end
 *	@param v_two the vertex at one end
 *	@param weight The weight of the edge
 */
void add_edge( vertex v_one, vertex v_two, edge_weight weight );

/**
 * @brief Runs the appropriate MST algorithm; each algorithm provides a
 * different implementation of this.
 *
 * @return the time it took the algorithm to run (in seconds)
 * @todo put the timing in the main program so that it is independent of the
 * algorithm and guaranteed to be uniform 
 */
double start_mst( void );

#endif /* __NEXTGEN_H */

/*  [Last modified: 2021 04 27 at 20:40:55 GMT] */
