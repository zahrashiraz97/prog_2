# prog_2
Program 2 of CSC 505

Programs must run from the command line and read input from standard input in gph format, as described here.
        c comment line 1
        ...
        c comment line k
        g number_of_nodes number_of_edges
	  --------------- (not part of the input)
        OPTIONAL node position info (for conversion to graphml)
        n v_1 x_1 y_1
        ...
        n v_n x_n y_n
        ------------
        ALWAYS
        e source_1 target_1 weight_1
        ...
        e source_m target_m weight_m

    v_1 through v_n are node numbers, typically 1 through n
    x_i, y_i are x and y coordinates of v_i


Sample program: 

Compilation: g++ -o file_name prims.cpp
Execution: graph.gph | ./file_name

Or:

./file_name
g 7 6
e 1 2 8
e 1 3 7
e 2 4 4
e 2 5 3
e 3 6 2
e 3 7 1
