# prog_2
Program 2 of CSC 505

We have implemented one version of Prim's and two versions of Kruskal's algorithms:
1. Prim's algorithm with lazy deletion (uses min heap)
2. Kruskal's with use of min heap and disjoint sets
3. Kruskal's where the edges are sorted 

Sample program: 

Compilation: g++ -o file_name program.cpp

Execution: cat graph.gph | ./file_name

Or:

./file_name
g 7 6

e 1 2 8

e 1 3 7

e 2 4 4

e 2 5 3

e 3 6 2

e 3 7 1
