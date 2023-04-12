/**
 *	@file	random_graph.c
 *	@brief	Creates a randomly generated, connected graph
 *	@author	Barry Peddycord
 *	@date	November 25, 2009
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
 * $Id: random_graph.c 273 2010-08-11 21:26:03Z mfms $
 */

/**
 * @todo make use of functions in graph module so that all graph
 * generation is consistent
 *
 * @todo some problems with cal calculation of expected collisions, gross
 * underestimate for dense graphs, overestimate for sparse ones.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>

#include <limits.h>
#include <math.h>
#include <time.h>

#include "common.h"
#include "nextgen.h"
#include "random.h"
#include "time.h"

/**
 * maximum edge weight, to be calculated based on number of vertices
 * to guarantee no overflow for a minimum spanning tree or path
 */
static edge_weight max_weight;

/**
 * maximumum number of collisions; to prevent thrashing in case of very dense
 * graphs; could be based on number of vertices, but ...
 */
#define COLLISION_LIMIT ((uint) 1 << 31)

/**
 * Hash values used to evaluate polynomials in polynomial hashing. One
 * determines bit position within a byte, the other determimes the byte.
 *
 * @todo There appears to be an odd interaction when asking for 25 vertices
 * and 129 edges with a seed of 1. It takes forever (may not terminate). But
 * other seeds appear to have a problem as well. The 25,128 combination is
 * okay. A different bit hash does not help, but a different byte hash value
 * does!
 */

#define BIT_HASH_VALUE 37
#define BYTE_HASH_VALUE 113

/**
 * Load factor for the bit vector, roughly the likelihood of a collision
 * between distinct edges; also determines the size of the bit vector; small
 * load factor implies large bit vector.
 *
 * Load factor determines the number of bytes, which will be number of edges
 * divided by load factor.
 *
 * In addition, within each byte, a potential edge is hashed to one of
 * HASH_BITS_IN_BYTE of its bits. 
 */
#define LOAD_FACTOR 0.05
#define HASH_BITS_IN_BYTE 7

/** 
 * The following ensure that there's some terminal output periodically
 */
#define TICK_ITERATIONS 1000000
static int ticker = TICK_ITERATIONS;

static void tick_if_needed( void ) {
  if ( ticker-- == 0 )
    {
      fprintf( stderr, ".");
      ticker = TICK_ITERATIONS;
    }
}

/**
 * for vertices
 */
#define SWAP( v_1, v_2 ) { uint tmp = v_1; v_1 = v_2; v_2 = tmp; }

/**
 * A utility function for computing the expected number of collisions.
 * @param b a number of balls to be tossed randomly into bins
 * @param B the number of bins, double because this number could be very
 * large, as in n (n - 1 ) / 2; b is double, too, for compatability
 * @param F the number of bins that already have a ball in them (also double
 * for compatibility)
 * @return the expected number of times a ball is tossed into a bin that
 * already has one or more balls in it; ech ball is tossed untl it lands in
 * an empty bin
 *
 * Let p_i be the probability that the i-th ball lands in a non-empty bin,
 * which is
 *      1 - (B - F) / B * ((B - 1) / B)^{i-1}
 * - in order to land in an empty bin, that bin has to be one of the B-F that
 * did not have a ball initially, and none of the previous i-1 balls landed
 * there.
 * The expected number of times we have to toss the ball before it lands in
 * an empty bin is
 *      p_i + (p_i)^2 + (p_i)^3 + ...
 *      = p_i / (1 - p_i)
 * or, in this case,
 *      B / (B - F) * (B / (B - 1))^{i-1} - 1
 * then the desired value is
 *      sum[ i = 1 to b ] p_i / (1 - p_i)
 * which becomes
 *      (B / (B - F)) * (B - 1) * ((B / (B - 1))^b - 1) - b
*/
static double balls_in_bins_collisions( double b, double B, double F ) {
#ifdef BDEBUG
  printf( "-> balls_in_bins_collisions: b = %f, B = %f, F = %f\n",
          b, B, F );
#endif
  
  assert( B - F > 0 );
  double base = B / (B - 1);
  double power = pow( base, b );
  double inner_sum = 1 - B - power + B * power;
  double result = B / (B - F) * inner_sum - b;
#ifdef BDEBUG
  printf( "  base = %f, power = %f, inner_sum = %f, result = %f\n",
          base, power, inner_sum, result );
#endif
  return result;
}

/**
 * Adds an edge if its hash value does not collide with a previously added
 * edge; keeps track of added edges in the edge_bit_vector. Assumes the edge
 * is not a self-loop.
 *
 * @param v_one one of the endpoints
 * @param v_two the other endpoint
 * @param weight the weight of the edge
 * @param force true if the edge should be added regardless of hash value
 * @return true if a new edge was generated/written
 *
 * Note: Always puts the lower numbered vertex first in case this needs to be
 * used to generate graphs for crossing minimization.
 *
 * Additional note: the hash value might have false positives, i.e., there
 * may be collisions in cases where two proposed edges are not identical.
 */
bool write_edge( vertex v_one, vertex v_two, edge_weight weight, bool force );

/**
 * Creates an initial random tree that includes all the vertices -- to ensure
 * that the graph will be connected.
 */
void create_initial_tree( vertex num_vertices );

/**
 * Uses a polynomial hash function to compute a hash value for a sequence of
 * bytes. If the k bytes are b_{k-1}, ..., b_0 the polynomial P(x) is
 * (b_{k-1} x^{k-1} + ... + b_0 x^0) mod z. The hash function evaluates P(x)
 * at x = some suitably chosen prime number and maps the result to an integer
 * in the range [0,max-1]; the modulus z should be >= max and relatively
 * prime to x (otherwise some integers in the range will never be mapped to)
 *
 * @param start_ptr pointer to b_{k-1}, the first byte in the sequence
 * @param length number of bytes in the sequence, i.e., k
 * @param x the prime number at which P(x) is evaluated
 * @param range_size the number of integers in the range, i.e., max
 * @return an (somewhat random) integer in the interval [0,max-1]; a
 * particular sequence of bytes always gives the same return value
 */
unsigned long polynomial_hash( const char * start_ptr,
                               unsigned int length,
                               unsigned int x,
                               unsigned long range_size );

/**
 * This is just a bit vector representing the *hash values* of the edges that
 * have already been added. Only one of two distinct edges is added if they
 * have the same hash value. The bit vector obviousky cannot be long enough
 * to accomodate every possible pair of vertices.
 */ 
static unsigned char * edge_bit_vector;

/**
 * @brief size (number of bytes) of the bit vector keeping track of hash
 * values of edges - edge_bit_vector */
static size_t bit_vector_length;

/**
 * @brief Constructs a random *connected* graph: edges are chosen randomly
 * and given uniformly random weights in the range [0,max_weight-1];
 * connectivity is ensured by constructing an initial spanning tree.
 *
 * @param num_verts number of vertices
 * @param num_edges number of edges
 */
void random_graph( size_t num_verts, size_t num_edges ) {
    /* guarantee that no MST nor shortest path will have total weight
       exceeding LONG_MAX */
    max_weight = (edge_weight) (MAX_EDGE_WEIGHT / (double) num_verts) + 1;
  uint collisions = 0;     /*< Number of collisions during hashing  */

  // ensure that the graph is at least a tree
  size_t edges_left = (num_edges > num_verts-1) ? num_edges : num_verts - 1;
	
  fprintf( stderr, "Building a graph with %zu vertices and %zu edges.\n",
           num_verts, edges_left );
  
#ifndef DIRECT_INPUT     /* write graph to to file */
	/* Open and write header data to the output file. */
    if ( ofname != NULL )
      outfile = fopen( ofname, "w" );
    else
      outfile = stdout;
	fprintf( outfile, "c random_graph: %zu vertices, %zu edges, seed = %d",
             num_verts, edges_left, seed );
  time_t raw_time = time(NULL);
  struct tm * current_time = gmtime(&raw_time);
  fprintf( outfile, ", %s", asctime(current_time));
	fprintf( outfile, "g %zu %zu\n", num_verts, edges_left );
#else  /* generate an empty graph with appropriate preamble internally */
    init_graph( num_verts, edges_left );
#endif

#ifdef DEBUG
    printf( "start initialization, num_verts = %zu,"
            " num_edges = %zu, edges_left = %zu\n",
            num_verts, num_edges, edges_left );
#endif

    size_t bytes_per_edge = (size_t) round( 1 / LOAD_FACTOR );
    // this should really be a prime number, but ...
    bit_vector_length = edges_left * bytes_per_edge + 1;
	edge_bit_vector = calloc( bytes_per_edge, edges_left );
    assert( edge_bit_vector != NULL );

  create_initial_tree( num_verts );
  edges_left -= num_verts - 1;

  double d_verts = (double) num_verts;
  double d_edges = (double) edges_left;
  double d_bits = (double) bit_vector_length * HASH_BITS_IN_BYTE;

  double expected_duplicates
    = balls_in_bins_collisions( d_edges, d_verts * (d_verts - 1) / 2,
                                d_verts - 1 );

  double expected_false_positives
    = balls_in_bins_collisions( d_edges, d_bits, d_verts - 1 );

  double expected_collisions
    = expected_duplicates + expected_false_positives;

  fprintf( stderr, "\nLOAD_FACTOR = %f, bits = %10.0f\n"
           " expected_collisisons = %f,"
           " expected_duplicates = %f,\n"
           " expected_false_positives = %f\n",
           LOAD_FACTOR, d_bits,
           expected_collisions, expected_duplicates, expected_false_positives );

#ifdef DEBUG
    printf( "initialization completed, edges_left = %zu,"
            " char_size = %lu, bytes_per_edge = %lu, bits_per_edge = %lu\n",
            edges_left, bytes_per_edge, bits_per_edge );
#endif

    /*
      while there are still edges to be added, pick two random vertices and
      attempt to connect them with a random weight
     */
    while ( edges_left > 0 ) {
#ifdef DEBUG
      printf( "new edge? edges_left = %zu\n", edges_left );
#endif
      vertex v_one = (vertex) genrand_int32() % num_verts + 1;
      vertex v_two = (vertex) genrand_int32() % num_verts + 1;
      if ( v_one == v_two ) continue;
      edge_weight weight = (edge_weight) (genrand_int32() % max_weight);
      bool success = write_edge( v_one, v_two, weight, false );
      if ( success ) {
        edges_left--;
#ifdef DEBUG
        printf( " success, edges_left = %zu\n", edges_left );
#endif
      }
      else {
        collisions++;
#ifdef DEBUG
        printf( " failure, collisions = %u\n", collisions );
#endif
      }

      // check to make sure this is not a hopeless venture:
      if ( collisions > COLLISION_LIMIT ) {
        fprintf( stderr, "\nToo many collisions: %u\n", collisions );
        fprintf( stderr, " either density is too great"
                 " or need a different seed, seed = %d\n", seed );
        break;
      }
    }

    fprintf( stderr,
             "\nDone adding edges: collisions = %d, bit_vector_size = %6.1f MB"
             "\n",
             collisions,
             ((double) bit_vector_length) / 1000000.0 );
	
#ifndef DIRECT_INPUT
    fclose( outfile );
#endif
	
	free( edge_bit_vector );
}

/* 
   Add an edge if its hash value does not
   collide with a previously added edge;
 */
bool write_edge( vertex v_one, vertex v_two, edge_weight weight, bool force ) {
#ifdef DEBUG
  printf( "-> write_edge: %zu, %zu, %d, %d\n", v_one, v_two, weight, force );
#endif
  assert( v_one != v_two );
  tick_if_needed();

  /* one entry for each vertex to comprise the "string" to be hashed */
  /* set up hash string so that the larger vertex number comes first */
  uint hashstr[2];
  hashstr[0] = (v_one > v_two) ? v_one : v_two;
  hashstr[1] = (v_one < v_two) ? v_one : v_two;
  /* get the bit position and the byte using polynomial hashing*/
  unsigned char bit_position
    = polynomial_hash( (char *) hashstr, 2 * sizeof(uint), BIT_HASH_VALUE,
                       HASH_BITS_IN_BYTE );
  uint byte_position
    = polynomial_hash( (char *) hashstr, 2 * sizeof(uint),
                       BYTE_HASH_VALUE, bit_vector_length );

  /* check for collision (and forcing) */
  unsigned char bit = 1 << bit_position; 
  if ( ( edge_bit_vector[ byte_position ] & bit )
       && ! force ) {
#ifdef DEBUG
      printf( "<- write_edge: collision, byte = %d, bit = %d\n",
              byte_position, bit_position );
#endif
      return false;
  }
  /* Mark the bit in the bit vector and add the edge */
  edge_bit_vector[ byte_position ] |= bit;
  actual_number_of_edges++;

  if ( v_one > v_two ) SWAP( v_one, v_two )

#ifndef DIRECT_INPUT
  fprintf( outfile, "e %zu %zu "EDGE_WEIGHT_FORMAT"\n", v_one, v_two, weight );
#else
  add_edge( v_one, v_two, weight );
#endif

#ifdef DEBUG
  printf( "<- write_edge: success\n" );
#endif
  return true;
}

void create_initial_tree( size_t num_vertices ) {
  // attach each vertex in turn to a random predecessor and choose a random
  // weight for the corresponding edge.
  for ( vertex v = 2; v <= num_vertices; v++ ) {
    vertex w = (genrand_int32() % (v - 1)) + 1;
    edge_weight weight = (edge_weight) (genrand_int31() % max_weight);
    write_edge( v, w, weight, true );
  }
}

/*
 * Uses a polynomial hash function to compute a hash value for a sequence of
 * bytes. If the k bytes are b_{k-1}, ..., b_0 the polynomial P(x) is
 * (b_{k-1} x^{k-1} + ... + b_0 x^0) mod z. The hash function evaluates P(x)
 * at x = some suitably chosen prime number and maps the result to an integer
 * in the range [0,max-1]; the modulus z should be >= max and relatively
 * prime to x (otherwise some integers in the range will never be mapped to)
 *
 * @param start_ptr pointer to b_{k-1}, the first byte in the sequence
 * @param length number of bytes in the sequence, i.e., k
 * @param x the prime number at which P(x) is evaluated
 * @param range_size the number of integers in the range, i.e., max
 * @return an (somewhat random) integer in the interval [0,max-1]; a
 * particular sequence of bytes always gives the same return value
 */
unsigned long polynomial_hash( const char * start_ptr,
                               unsigned int length,
                               unsigned int x,
                               unsigned long range_size )
{
#ifdef DEBUG
  printf( "-> polynomial hash: length = %u, x = %u, range_size = %lu\n",
          length, x, range_size );
#endif
  unsigned int modulus = range_size;
  /* adjust the modulus to ensure that it's relatively prime to evaluate at */
  if ( modulus % x == 0 || x % modulus == 0 )
    modulus++;
  /* evaluate the polynomial using Horner's rule */
  const char * current_byte = start_ptr;
  unsigned long hash_value = 0;
  while ( current_byte < start_ptr + length ) {
    unsigned int coefficient = (unsigned int) (* current_byte);
    hash_value = (x * hash_value + coefficient) % modulus;
    current_byte++;
  }
#ifdef DEBUG
  printf( "<- polynomial hash: hash_value = %lu, range_size = %lu\n",
          hash_value, range_size );
#endif
  return hash_value % range_size;
}
