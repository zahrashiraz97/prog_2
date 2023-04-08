These are programs and scripts for generating
 four different types of graphs originally used in minimum spanning tree experiments.

To compile everything, simply say `make`

**Suggestion:** Use the `generate.sh` wrapper script to generate graphs of all the available types.
This guarantees a uniform approach. Each generated file has a standardized name
and the correct information about number of edges on the `g` line.
The first line always has the form
```
c NUM_NODES NUM_EDGES MAX_WEIGHT SEED
```
Programs can harvest this information if, for example, it turns out to be useful for creating csv compatible output.

random - completely random

geometric - vertices are random points in a unit square, edges are between vertices that are within distance d from each other, using the infinity norm: d(p_1,p_2) = max(|p_1.x - p_2.x|, |p_1.y - p_2.y|)

geo_wrap - like geometric but with wraparound, i.e., points near the bottom boundary are treated as close to those near the top boundary; and those near the left are close to those near the right

Delaunay triangulations - vertices are points in a unit square; edges form triangles; no point is inside the "circumcircle" of any triangle, i.e., the circle that includes all three triangle points

To generate these in .gph format, run, for the first three types,
   ./nextgen type n m seed file
where type is random, geometric, or geo_bound, n and m, are number of vertices and edges, respectively, seed is an integer random seed, and file is the output file. The result is always connected. With random graphs the number of edges is exactly the desired number; with geometric and geo_bound, it is close for very large graphs, but may be less than desired for smaller ones.

For the ***small*** triangulations, run
    ./random_delaunay.py
The -h option will give you details. Anything beyond 1000 nodes will take minutes and runtime appears to be quadratic.

For larger ones, use
```triangulation NUMBER_OF_NODES SEED [OUPUT_FILE]```
not as many options (infinite face of duals), but runs in linear time and can generate a triangulation with 10,000,000 nodes in about a minute.

Do 'make' to create the programs.

Use the gph2cnf script to convert the files generated by the above programs from their gph format to the cnf format required by some vertex cover programs. 