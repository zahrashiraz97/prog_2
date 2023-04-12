/**
 *	@file	nextgen.c
 *	@brief	Creates a randomly generated, connected graph and saves it as
 *			a .gph file.
 *	@author	Barry Peddycord
 *	@date	November 25, 2009
 * @modified 2010/06/11, Tony Chen: ensures connectivity for geometric graphs
 *	with some loss of "locality"
 *
 * [Part of Prof. Stallmann's NCSU MST project]
 *
 * This program accepts a number of vertices, a density value, and a random
 * seed, and then outputs a connected graph for use in MST searching. The
 * graphs generated are completely abstract with arbitrarily chosen edge
 * weights.
 *
 * The generator runs in linear time, running at best in linear time around
 * the number of edges in the graph, plus the time spent dealing with
 * collisions. The system uses a Hash function to test for collisions, with
 * a maximum of a 12% load factor.
 *
 * If compiled with DIRECT_INPUT, rather than writing the graph to a file,
 * it will push the data directly into a graph algorithm.
 *
 * Even though undirected for MST experiments, all graphs are dags if vertex
 * order for each edge is used: order is upper-left to lower-right for
 * geometric graphs and based on vertex number for random graphs.
 *
 * $Id: nextgen.c 262 2010-08-09 21:21:54Z mfms $
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>

#include <libgen.h>             /* basename() */
#include <limits.h>
#include <math.h>

#include "common.h"
#include "nextgen.h"
#include "random.h"
#include "time.h"

/** maximum length of a file name */
#define NAME_LENGTH 511

/** return value if too many edges */
#define TOO_MANY_EDGES 2

#ifndef DIRECT_INPUT

FILE *outfile;				/*< The output file for storing a random graph */

#else // DIRECT_INPUT

static double runningtime;			/*< The runtime of the algorithm. */
static FILE *reportfile;			/*< The file to save report data. */

#endif

uint seed;                /*< User-chosen random seed. */
char * ofname;            /*< The output filename. (NULL if stdout) */
size_t num_verts;         /*< Number of vertices. */
size_t num_edges;         /*< Number of edges. */
/** 
 * actual number of edges when the graph has been generated 
 * (may be different than asked for with a geometric graph)
 */
size_t actual_number_of_edges = 0;


/**
 * @brief Parses the command line options.
 *
 * @return true if the number of arguments is correct; prints a usage message
 * and returns false otherwise
 */
bool handle_args( int argc, char **argv )
{
  /* Checking for command line parameters. */
    if ( argc < 5 || argc > 6 ) {
        const char * from_where = argv[0];
        fprintf( stderr, "Usage: %s TYPE n m s OFILE\n", from_where );
        fprintf( stderr, " where TYPE = random, geometric or geo_wrap, n = # vertices, \n");
        fprintf( stderr, " m = # edges, s = seed, OFILE = output filename\n" );
        fprintf( stderr, " geo_wrap is geometric with wraparound (nodes close to a boundary connect with those on opposite side)\n" );
        return false;
    }
    else {
        sscanf( argv[2], "%zu", &num_verts );
        sscanf( argv[3], "%zu", &num_edges );
        sscanf( argv[4], "%u", &seed );
        if ( argc == 6 ) {
            ofname = argv[5];
        }
        else {
            ofname = NULL;
        }
    }
    return true;
}

/**
 * @brief This main program produces a random graph and serves two distinct
 * functions:
 *
 * -# if DIRECT_INPUT is defined the random graph is fed directly into an
      algortithm, the one whose implementation is linked with this module

   -# otherwise, the program writes the random graph to a file
 *
 * @todo separate these two functions into distinct programs or at least let
 * the user choose to output the randomly generated graph via a command-line
 * option.
 */
int main( int argc, char **argv ) {
    if ( ! handle_args( argc, argv ) )
      return EXIT_FAILURE;

    // abort if number of edges if more than C(n,2) and issue warning if
    // number of edges is too few to connect the graph
    unsigned long long int max_edges =
      (unsigned long long int) num_verts * (num_verts - 1) / 2;
    if ( num_edges >  max_edges ) {
        fprintf(stderr, "Too many edges requested, %zu; maximum possible number"
                " for %zu vertices = %llu\n",
                num_edges, num_verts, max_edges);
        return TOO_MANY_EDGES;
    }
    if ( num_edges < num_verts - 1 ) {
      printf( "Warning: Number of edges < number of vertices - 1\n" );
      printf( "Edges will be added to form a tree\n" );
    }

	init_genrand( seed );
	
    bool success = true;        /* geometric graphs may not be connected */
	if ( strcmp(argv[1], "random" ) == 0 ) {
      random_graph( num_verts, num_edges );
    }
    else if ( strcmp(argv[1], "geo_wrap" ) == 0 ) {
      wrap_around = true;
      success = geometric_graph( num_verts, num_edges );
    }
    else if ( strcmp(argv[1], "geometric" ) == 0 ) {
      wrap_around = false;
      success = geometric_graph( num_verts, num_edges );
    }
    else {
      fprintf( stderr, "Unknown graph type %s\n", argv[1] );
    }

    if ( ! success ) {
        fprintf(stderr,
                "*** WARNING: geometric graph is not connected:\n"
                "             vertices = %zu, desired edges = %zu ***\n",
              num_verts, num_edges);
        return EXIT_FAILURE;
    }

    fprintf(stderr, "+++ graph generation successful +++\n" );
#ifdef DIRECT_INPUT
    runningtime = start_mst();
    if ( ofname != NULL )
      reportfile = fopen( ofname, "a" );
    else reportfile = stdout;
    const char * base_name = basename( argv[0] );
    const char * graph_type = argv[1];
    char time_buffer[ BUFFER_SIZE ];
    strftime( time_buffer, sizeof( time_buffer ), "%F %T", time( NULL ) )
    fprintf( reportfile, "%s,%s,%u,%u,%lf,%llu,%s,%u\n",
             base_name, graph_type,
             num_verts, actual_number_of_edges,
             runningtime, COMPARISONS, time_buffer, seed );
    fclose( reportfile );
#endif
	return EXIT_SUCCESS;
}

/**
 * @todo Need a separate clock module for time-related issues (right now each
 * algorithm has its own).
 */

/* void timestamp( const char * message ) { */
/*   printf("[%7.3f]  %s", */
/*          (clock() - start_clock) / (double) CLOCKS_PER_SEC, message); */
/* } */
