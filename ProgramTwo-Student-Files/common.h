/*******************************************************************
 * Barry Peddycord
 *
 * 	Common data shared between all files.
 ******************************************************************/

#ifndef __COMMON_H
#define __COMMON_H

#include <stdlib.h>
#include <stdbool.h>

#define BUFFER_SIZE 512         /* holds standard Unix line - suffient for
                                   most purposes */

/**
 * Increment used when dynamically allocating an array that expands at runtime
 */
#define ALLOCATION_INCREMENT 16

typedef unsigned int uint;

/**
 * @todo Decide whether the number of vertices/edges should be global or
 * not. They have different meaning at the start generating random geometric,
 * where num_edges is just an estimate. 
 */
extern size_t num_verts;
extern size_t num_edges;

extern unsigned long long int comparisons;

#define INIT_COMPARISONS (comparisons = 0)
#define COMPARISONS comparisons

#define INCREMENT_COMPARISONS (comparisons++)

/**
 * for geometric graphs -- see geometric.c (should really be declared in a
 * separate geometric.h header, but none exists)
 */
extern bool wrap_around;

#endif /* __COMMON_H */

/*  [Last modified: 2021 04 27 at 20:37:58 GMT] */
