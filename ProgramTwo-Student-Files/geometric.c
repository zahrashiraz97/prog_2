/**
 *	@file	geometric.c
 *	@brief	Creates a randomly generated graph based on geometric considerations
 *	@author	 Matt Stallmann, based on work of Shengyen Tony Chen
 *	@date	July 2010
 *
 * [Part of Prof. Stallmann's NCSU MST project]
 *
 * Creates a random graph whose vertices are mapped to points in the plane
 * and whose edges are determined by their distance from other vertices. The
 * weight of each edge is its distance. Distances are calculated using the
 * 'infinity' or 'max' norm
 *
 * If compiled with DIRECT_INPUT, rather than writing the graph to a file,
 * it will push the data directly into a graph algorithm.
 *
 * @todo reexamine the use of DIRECT_INPUT and consider using graph.[ch]
 * @todo the mechanism for adding extra edges to guarantee
 * connectivity is arcane and does not work
 * @todo would be better to use Euclidian distance rather than the
 * infinity norm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include <stdbool.h>
#include <limits.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "common.h"
#include "nextgen.h"
#include "random.h"
#include "time.h"

/**
 * true if the graph is effectively embedded on a torus, i.e., points close
 * to the bottom edge are close to those near the top edge and ditto for left
 * and right edges.
 */
bool wrap_around = true;

/**
 * if there's no wrap around we need to account for vertices near the
 * boundaries when estimating the distance threshold that will get close to
 * the right number of edges (not clear that this is needed)
 */
#define BOUNDARY_FUDGE_FACTOR 0

/**
 * The following ensure that there's some terminal output periodically
 */
#define TICK_ITERATIONS 1000000
static int ticker = TICK_ITERATIONS;

static void tick_if_needed( void ) {
  if ( ticker-- == 0 )
    {
      fprintf(stderr, ".");
      ticker = TICK_ITERATIONS;
    }
}

#ifndef DIRECT_INPUT

/* #include "graph.h" */
/* static GRAPH * geo_graph; */

#endif

typedef struct point {
  int x;
  int y;
} POINT;

#define UP_AND_LEFT( p1, p2 ) ((p1.y < p2.y) || (p1.y == p2.y && p1.x < p2.x)) 

#define MIN( a, b ) (( a > b ) ? b : a) 
#define MAX( a, b ) (( a > b ) ? a : b)

/**
 * for vertices
 */
#define SWAP( v_1, v_2 ) { vertex tmp = v_1; v_1 = v_2; v_2 = tmp; }

/**
 * distance between two points p1, p2; use infinity norm and wrap-around
 */
static uint distance( POINT p1, POINT p2 ) {
#ifdef DEBUG
  printf( "-> distance: (%7.4f, %7.4f) , (%7.4f, %7.4f)\n",
          p1.x / (double) INT_MAX,
          p1.y / (double) INT_MAX,
          p2.x / (double) INT_MAX,
          p2.y / (double) INT_MAX );
#endif
  uint x_diff_normal = abs( p1.x - p2.x );
  uint x_diff_wrap = INT_MAX - x_diff_normal;
  uint x_diff = MIN( x_diff_normal, x_diff_wrap );
  uint y_diff_normal = abs( p1.y - p2.y );
  uint y_diff_wrap = INT_MAX - y_diff_normal;
  uint y_diff = MIN( y_diff_normal, y_diff_wrap );
#ifdef DEBUG
  printf( "  distance: xdiff = %7.4f, ydiff = %7.4f\n",
          x_diff / (double) INT_MAX,
          y_diff / (double) INT_MAX );
#endif
  return MAX(x_diff, y_diff);
}

/**
 * @return the Euclidean distance between the points, to be used to
 * compute edge weights
 */
static double euclidian_distance(POINT p1, POINT p2) {
    double xdiff = p1.x - p2.x;
    double ydiff = p1.y - p2.y;
    return sqrt(xdiff * xdiff + ydiff * ydiff);
}

/** Creates n disjoint sets numbered 1, ... ,n */
static void make_sets( uint n );

/** returns the canonical element (number) of the set to which v belongs */
static vertex find_set( vertex v );

/**
 * merges the two sets to which v1 and v2 belong, if necessary.
 * @return true if the sets were disjoint (v1 and v2 were in different sets)
 * initially
 */
static bool set_union( vertex v1, vertex v2 );

/**
 * @return the current number of disjoint sets
 */
static uint get_number_of_sets( void );

/** deallocates structures associated with disjoint sets */
static void delete_sets( void );

/* typedef enum direction_enum { */
/*   DOWN, */
/*   RIGHT, */
/*   DOWN_LEFT, */
/*   DOWN_RIGHT, */
/*   NUMBER_OF_DIRECTIONS */
/* } DIRECTION; */

typedef struct square {
  vertex * vertex_list;
  uint number_of_vertices;
  /* keep track of closest unconnected pair in case extra edges are needed */
  vertex closest_unconnected_v1;
  vertex closest_unconnected_v2;
  uint closest_unconnected_distance;
} SQUARE;

static POINT * vertex_coordinates = NULL;
static vertex largest_vertex_number = 0;

/**
 * @return true if the coordinates of v_1 are above and/or (if the same
 * y-coordinate) to the left of those of v_2
 */
static bool up_and_left( vertex v_1, vertex v_2 ) {
  return UP_AND_LEFT( vertex_coordinates[ v_1 ], vertex_coordinates[ v_2 ] );
}

/**
 * The grid is logically stored as a two-dimensional array, but physically as
 * a one-dimensional array of SQUARES; addressing in the one dimensional
 * array is standard: [x][y] -> grid_dimension * x + y
 */
SQUARE * grid_squares = NULL;
/** number of squares in each row or column */
uint grid_dimension = 0;
/** total number of squares = grid_dimension * grid_dimension */
uint number_of_squares = 0;

/** two vertices whose points are within radius of each other will be
    connected */
uint radius = 0;
/** radius expressed as a "real" number; currenly used for more accuracy in
    grid square calculations */
double real_radius = 0.0;

/**
 * Initializes all global values and date structures
 * @param n number of vertices (also largest vertex number)
 */
static void init_geometric_graph( vertex n );

/**
 * Deallocates all structures
 */
static void delete_geometric_graph( void );

/**
 * Chooses a random non-negative integer point for each vertex and adds the
 * vertex to the appropriate grid square. If output is going to a file, the
 * coordinates of each point are output in format
 *     n id x y
 */
static void put_vertices( void );

/**
 * Adds all the edges
 */
static void add_all_edges( void );

/**
 * Adds all of the edges among vertices in the same grid square, identified
 * by its index.
 */
static void add_edges_in_square( uint square_index );

/**
 * Adds all of the edges between vertices in two distinct (adjacent) grid
 * squares.
 * @param square_one_index index of the square *from* which edges are to be
 * added (this is upward and/or to the left).
 * @param square_two_index index of the square *from* which edges are to be
 * added (this is downward and/or to the right).
 */
static void add_edges_between_squares( uint square_one_index,
                                       uint square_two_index );

/**
 * Adds extra edges to the graph to ensure connectivity; this is done by
 * choosing the lowest cost unconnect pair from each square; squares are
 * permuted randomly before this is done.
 */
static void add_extra_edges( void );

/**
 * Randomly permutes the square in the grid; Caution: the grid is no longer
 * valid after this function is called.
 */
static void permute_squares( void );

/**
 * Creates a geometric graph - edges are added via add_edge() if the result
 * is to be fed directly into an algorithm (DIRECT_INPUT) or printed to a
 * file if not.
 * @return true if the graph is connected
 */
bool geometric_graph( size_t num_verts, size_t num_edges ) {

  init_geometric_graph( num_verts ); /* sets largest_vertex_number */

#ifndef DIRECT_INPUT
  // geo_graph = new_graph( largest_vertex_number, num_edges );
  // add comments (header for graph output file will come later since we
  // don't know number of edges yet).
  if ( ofname != NULL ) {
    outfile = fopen( ofname, "w" );
    if ( outfile == NULL ) {
      printf( "Unable to open file %s for writing\n", ofname );
      exit(1);
    }
  }
  else outfile = stdout;
  char * type_string = "geometric";
  char * wrap_around_string = "";
  if ( wrap_around ) {
    type_string = "geo_wrap";
    wrap_around_string = " with wraparound";
  }
  fprintf( outfile, "c random geometric graph%s: %zu vertices, %zu (desired) edges,"
             " seed = %d\n",
           wrap_around_string, num_verts, num_edges, seed );
  fprintf(outfile, "c created by 'nextgen %s %zu %zu %d'",
          type_string, num_verts, num_edges, seed );
  time_t raw_time = time(NULL);
  struct tm * current_time = gmtime(&raw_time);
  fprintf( outfile, ", %s", asctime(current_time));
  fprintf( outfile, "g %zu %zu\n", largest_vertex_number, num_edges );
  /**
   * @todo create a mechanism for adding comments to a graph
   */
#else
  init_graph( largest_vertex_number, num_edges );
#endif
  put_vertices();
  fprintf( stderr, "Done putting vertices\n" );
  make_sets( num_verts );
  add_all_edges();
  size_t components = get_number_of_sets();
  fprintf( stderr, "\n");
  fprintf( stderr, "- Number of regular edges = %zu\n", actual_number_of_edges );
  fprintf( stderr, "- Number of components before extra edges = %zu\n", components );

  add_extra_edges();

  components = get_number_of_sets();
  fprintf(stderr, "\n");
  fprintf( stderr, "- Number of edges after extra edges = %zu\n",
          actual_number_of_edges );
  fprintf( stderr, "- Number of components after extra edges  = %zu\n", components );

  fprintf( stderr, "Done adding edges\n" );
  delete_sets();
  delete_geometric_graph();
#ifndef DIRECT_INPUT

/*   outfile = fopen( ofname, "w" ); */
/*   if ( outfile == NULL ) { */
/*     printf( "Unable to open file %s for writing\n", ofname ); */
/*     exit(1); */
/*   } */
/*   write_graph( geo_graph, outfile ); */
  fclose( outfile );
#endif
  if ( components > 1 ) return false;
  else return true;
}

static void init_geometric_graph( vertex n ) {
  largest_vertex_number = n;
  /* Allocate space for vertex information */
  vertex_coordinates
    = calloc( largest_vertex_number + 1, sizeof(POINT) );
  assert( vertex_coordinates != NULL );

  /*
    Use desired number of edges to calculate radius and number of radius x
    radius grid squares. Allocate memory for grid squares.

    let n_d = expected number of neighbors of a vertex within radius d;

    let m and n be the number of edges and vertices, respectively;

    then n_d = (2*d)^2 and the expected number of edges is (2*d)^2 * n/2 (an
    edge corresponds to two neighbors) so if we want m = (2*d)^2 * n/2, then
    we need d = sqrt( m / 2*n ) = sqrt( 2*m ) / 2*n
  */
  double ratio = sqrt( 2.0 * num_edges ) / (2.0 * largest_vertex_number);
  if ( ! wrap_around ) {
    ratio = ratio * (1 + BOUNDARY_FUDGE_FACTOR * ratio);
  }
  real_radius = INT_MAX * ratio;
  radius = (uint) real_radius;
  grid_dimension
    = ((uint) (INT_MAX / real_radius)) + 1; /* +1 needed due to truncation */
  
  number_of_squares = grid_dimension * grid_dimension;
  grid_squares
    = calloc( number_of_squares, sizeof(SQUARE) );
  assert( grid_squares != NULL );

  fprintf( stderr, "geo graph: ratio = %f, radius = %d, grid_dimension = %u\n",
         ratio, radius, grid_dimension );
 
  fprintf( stderr, "num_squares = %u\n", number_of_squares);

  // initialize grid square data so that closest unconnected distance is
  // "infinity"; and the direction for choosing closest unconnected points is
  // chosen randomly
  for ( uint i = 0; i < number_of_squares; i++ ) {
    grid_squares[i].vertex_list = NULL; /* not really needed, but ... */
    grid_squares[i].closest_unconnected_distance = INT_MAX;
/*     grid_squares[i].closest_unconnected_direction */
/*       = genrand_int32() % NUMBER_OF_DIRECTIONS; */
  }
}

static void delete_geometric_graph( void ) {
  free( vertex_coordinates );
  for ( uint i = 0; i < number_of_squares; i++ )
    free( grid_squares[i].vertex_list );
  free( grid_squares );
}

static void put_vertices( void ) {
#if defined(DEBUG)
  printf("\nVertex\t\tx\t\ty\tgrid square/square number\n");
#endif
  for( vertex i = 1; i <= largest_vertex_number; i++) {
    tick_if_needed();
    vertex_coordinates[i].x = genrand_int32() % INT_MAX;
    vertex_coordinates[i].y = genrand_int32() % INT_MAX;

#ifndef DIRECT_INPUT
    fprintf( outfile, "n "VERTEX_FORMAT" %d %d\n",
             i, vertex_coordinates[i].x, vertex_coordinates[i].y );
#endif

    /* Calculates the position of the grid square for the vertex. */
    uint square_x = (uint) (vertex_coordinates[i].x / radius);
    uint square_y = (uint) (vertex_coordinates[i].y / radius);
    uint index = square_x * grid_dimension + square_y;
    uint n_verts = grid_squares[ index ].number_of_vertices;
    if ( n_verts % ALLOCATION_INCREMENT == 0 ) {
      grid_squares[index].vertex_list
        // bug fix: added parens (mfms, 2010/07/21)
        = realloc( grid_squares[index].vertex_list,
                   (n_verts + ALLOCATION_INCREMENT) * sizeof(vertex) );
      assert( grid_squares[index].vertex_list != NULL );
    }
    grid_squares[index].vertex_list[n_verts] = i;
    grid_squares[index].number_of_vertices++;

#if defined(DEBUG)
    printf( "%u\t%10.4f\t%10.4f\t(%u,%u)\t%u\n", i,
            vertex_coordinates[i].x / (double) INT_MAX,
            vertex_coordinates[i].y / (double) INT_MAX,
            square_x, square_y, index );
#endif
  }
}

/**
 * Adds an edge to the graph, taking into account the different ways the
 * graph will be used (directly fed into an algorithm or printed to a file).
 * Also updates the disjoint sets appropriately.
 */
static void add_an_edge( vertex v1, vertex v2, edge_weight weight ) {
  tick_if_needed();
  if ( ! up_and_left( v1, v2 ) ) SWAP( v1, v2 );
#ifdef DIRECT_INPUT
  add_edge( v1, v2, weight );
#else
  fprintf(outfile, "e "VERTEX_FORMAT" "VERTEX_FORMAT" "EDGE_WEIGHT_FORMAT"\n",
          v1, v2, weight);
#endif
  actual_number_of_edges++;
  set_union( v1, v2 );
#if defined(DEBUG)
  printf( "add_edge( %u, %u, %7.4f )\n", v1, v2,
          weight / (double) INT_MAX );
#endif
}

static void add_all_edges( void ) {
  // for each grid square:
  //
  // connect all vertices within the square
  //
  // connect vertices within the square to those on the row below and
  // those to the right on the same, if they are close enough; wrap around if
  // needed;
  // in other words, from square (i,j), check squares (i+1,j-1), (i+1,j),
  // (i+1,j+1), and (i,j+1) 
  for ( uint square_index = 0;
        square_index < number_of_squares;
        square_index++ )
    {
      uint square_i = square_index / grid_dimension;
      uint square_j = square_index % grid_dimension;
      // within square
#if defined(DEBUG)
      printf( " -- connecting within: square = %u, i = %u, j = %u\n",
              square_index, square_i, square_j );
#endif
      add_edges_in_square( square_index );
      
      // determine direction for keeping track of closest unconnected point
/*       DIRECTION direction = grid_squares[ square_index ] */
/*         .closest_unconnected_direction; */

      // calculate 2-D indexes of four neighboring squares
      uint down_i = (square_i + 1) % grid_dimension;
      uint left_j
        = (square_j + grid_dimension - 1) /* must be positive */
        % grid_dimension;
      uint right_j = (square_j + 1) % grid_dimension;
      // The neighbors we want are [i,right_j]; [i+1,left_j]; [i+1,j]; and
      // [i+1,righ_j]

      // Caluclate indexes of the appropriate squares
      uint current = square_i * grid_dimension + square_j;
      uint right = square_i * grid_dimension + right_j;
      uint down_left
        = down_i * grid_dimension + left_j;
      uint down = down_i * grid_dimension + square_j;
      uint down_right
        = down_i * grid_dimension + right_j;
      // add the edges
      if ( wrap_around || square_j < grid_dimension - 1 ) {
        add_edges_between_squares( current, right );
#if defined(DEBUG)
      printf( " <> Connecting between squares [%u,%u] and [%u,%u]\n",
              square_i, square_j,
              square_i, right_j );
#endif
      }
      if ( wrap_around || ( square_i < grid_dimension - 1 && square_j > 0 ) ) {
        add_edges_between_squares( current, down_left );
#if defined(DEBUG)
        printf( " <> Connecting between squares [%u,%u] and [%u,%u]\n",
                square_i, square_j,
                down_i, left_j );
#endif
      }
      if ( wrap_around || square_i < grid_dimension - 1 ) {
        add_edges_between_squares( current, down );
#if defined(DEBUG)
        printf( " <> Connecting between squares [%u,%u] and [%u,%u]\n",
                square_i, square_j,
                down_i, square_j );
#endif
      }
      if ( wrap_around || ( square_i < grid_dimension - 1
                            && square_j < grid_dimension - 1 ) ) {
        add_edges_between_squares( current, down_right );
#if defined(DEBUG)
        printf( " <> Connecting between squares [%u,%u] and [%u,%u]\n",
                square_i, square_j,
                down_i, right_j );
#endif
      }
    }
}

static void add_edges_in_square( uint square_index ) {
  // for every pair of vertices in the square add an edge; since there are no
  // modifications here, we can use a copy of the square
  SQUARE square = grid_squares[ square_index ];
  for ( uint i = 0; i < square.number_of_vertices; i++ ) {
    for ( uint j = i + 1; j < square.number_of_vertices; j++ )
      {
        uint v_i = square.vertex_list[i];
        uint v_j = square.vertex_list[j];
        POINT p_i = vertex_coordinates[v_i];
        POINT p_j = vertex_coordinates[v_j];
        edge_weight weight = (edge_weight) round(euclidian_distance(p_i, p_j));
        add_an_edge( v_i, v_j, weight );
      }
  }
}

static void add_edges_between_squares( uint square_one_index,
                                       uint square_two_index ) {
  // the code for determining whether an edges should be created does not
  // require modification of square data; a copy of the squares are used here
  SQUARE square_one = grid_squares[ square_one_index ];
  SQUARE square_two = grid_squares[ square_two_index ];
  for ( uint i = 0; i < square_one.number_of_vertices; i++ ) {
    for ( uint j = 0; j < square_two.number_of_vertices; j++ )
      {
        uint v_i = square_one.vertex_list[i];
        uint v_j = square_two.vertex_list[j];
        POINT p_i = vertex_coordinates[v_i];
        POINT p_j = vertex_coordinates[v_j];
        tick_if_needed();
        uint inf_norm = distance(p_i, p_j);
        if ( inf_norm <= radius ) {
            edge_weight weight = (edge_weight) round(euclidian_distance(p_i, p_j));
            add_an_edge(v_i, v_j, weight);
        } // end, if weight <= radius
        else if ( //is_unconnected_direction &&
                  inf_norm < 
                  grid_squares[square_one_index]
                  .closest_unconnected_distance ) {
          // here we need to modify data in the actual square
          grid_squares[square_one_index].closest_unconnected_distance
            = inf_norm;
          grid_squares[square_one_index].closest_unconnected_v1
            = v_i;
          grid_squares[square_one_index].closest_unconnected_v2
            = v_j;
        } // end, weight < that of closest unconnected point
      }
  }
}

/**
 * traverses grid squares in random order until everthing is connected using
 * closest unconnected edges
 */
static void add_extra_edges( void ) {
    permute_squares();
    fprintf( stderr, "Done permuting squares\n");
    for ( uint i = 0; i < number_of_squares; i++ ) {
        if ( get_number_of_sets() == 1 ) break;
        if ( grid_squares[i].closest_unconnected_distance == INT_MAX ) continue;
        vertex v1 = grid_squares[i].closest_unconnected_v1;
        vertex v2 = grid_squares[i].closest_unconnected_v2;
        POINT p1 = vertex_coordinates[v1];
        POINT p2 = vertex_coordinates[v2];
        edge_weight weight = (edge_weight) round(euclidian_distance(p1, p2));
        add_an_edge(v1, v2, weight);
    }
}

static void permute_squares( void ) {
  SQUARE temp;
  uint i = number_of_squares;
  while ( --i > 0 ) {
    tick_if_needed();
    /* swap grid_squares[i] with a random element among grid_squares[0...i] */
    uint j = genrand_int32() % (i + 1);
    if ( j != i ) {
      temp = grid_squares[i];
      grid_squares[i] = grid_squares[j];
      grid_squares[j] = temp;
    }
  }
}


/*
  Implementation of disjoint sets using merge by rank and path compression
 */

typedef struct uf_set {
    vertex parent;
    vertex rank;
} UNION_FIND_SET;

static UNION_FIND_SET * ufset;
static uint number_of_sets;

/** Creates n disjoint sets numbered 1, ... ,n */
static void make_sets( uint n ) {
  number_of_sets = n;
  ufset = (UNION_FIND_SET * ) calloc( n + 1, sizeof(UNION_FIND_SET) );
  assert( ufset != NULL );
    for ( vertex i = 1; i <= n; i++ )
    {
        ufset[i].rank = 0;
        ufset[i].parent = i;
    }
}

/** returns the canonical element (number) of the set to which v belongs */
static vertex find_set( vertex v ) {
  if ( ufset[v].parent != v )
    {
      ufset[v].parent = find_set( ufset[v].parent );
    }
  return ufset[v].parent;
}

/**
 * merges the two sets to which v1 and v2 belong, if necessary.
 * @return true if the sets were disjoint (v1 and v2 were in different sets)
 * initially
 */
static bool set_union( vertex v1, vertex v2 ) {
  vertex a = find_set(v1);
  vertex b = find_set(v2);
  if ( a == b ) return false;
  if ( ufset[a].rank > ufset[b].rank )
    {
      ufset[b].parent = a;
    }
    else
    {
      ufset[a].parent = b;
      if ( ufset[a].rank == ufset[b].rank ) ufset[b].rank++;
    }
  number_of_sets--;
  return true;
}

/**
 * @return the current number of disjoint sets
 */
static uint get_number_of_sets( void ) { return number_of_sets; }

/** deallocates structures associated with disjoint sets */
static void delete_sets( void ) {
  free( ufset );
}

/*  [Last modified: 2021 04 29 at 19:09:26 GMT] */
