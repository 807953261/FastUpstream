																			Compiling

A Makefile is provided (the executable is giscup.out).

Dependency: a folder "libs" with the RapidJSON libraryis required.

To compile, just type "make" (the executable will be giscup.out)

Explanation:
The compiler flag: “-I libs/”  adds the RapidJson library
The g++ flag: “-fopenmp -D_GLIBCXX_PARALLEL”  enables OpenMP and the parallel mode of the STL library.
Furthermore, we used the -O3 optimization flag and  -DNDEBUG to disable the assertions/debugging mode.

------------------------------------------------------------------------------------------------------
																		Libraries used

We only used C++11 and RapidJson library (according to the GISCUP website, it is ok to use libraries 
for assisting with reading the input files). Furthermore, we used OpenMP (according to one of the 
messages in the GISCUP discussion list, g++ is available and OpenMP is allowed).

------------------------------------------------------------------------------------------------------
																				Running

The program expects three parameters (as described in the GISCUP website):
1) The path to the input JSON file
2) The path to the file with the starting points
3) The path to the output file

OBSERVATION 1: Since our implementation is I/O bound, it is much faster if the input file is already
in the disk caches (i.e., if it is executed twice the second execution will be much faster).

OBSERVATION 2: Since, according to the GISCUP organizers, the output does not have to be sort and may
contain repeated IDs we do not sort the results

------------------------------------------------------------------------------------------------------
														Main ideas behind the algorithm

- Our main algorithm is very simple (though the parser has A LOT of optimizations).
- We deal with a graph (instead of with the geometry) created based on the JSON file.
- Our algorithm is based on the idea of Biconnected Components and Block-Cut trees (BC-Trees).
- We use Tarjan’s algorithm to find the articulation points of the graph and, then, construct a 
block-cut tree.
- The result is found by quickly (in O(N) time) processing the block-cut tree (we basically have to 
analyze which block of the tree have “special” vertices (starting points/controllers)).
- To simplify the treatment of the special cases, our “solver” considers only vertices can be starting 
points. If we have an edge that is a starting point → we create an “artificial” vertex:
-- If edge e=(a,b) is an starting point --> we create an "artificial/dummy" vertex v and pretend v is 
a starting point.
-- Then, e=(a,b) is replaced with e'=(a,v) and e''=(v,b)   (if a starting point line is composed of 
several edges --> we do that for all edges)
- This is done during the creation of the graph (in the parsing step) --> thus, our "solver" does 
not have to handle this special case directly.

- Some notes (and figures) describing our idea can be found here: 
https://docs.google.com/document/d/1mGpnGiNkoCkDR7OVXr1feKB3KIAsGVNo8QLRRbvte0s/edit?usp=sharing

------------------------------------------------------------------------------------------------------
																			More details
- For performance, we employed several techniques including:
-- Usage of a fast custom parser.
-- Usage of parallel programming.
-- Overlapping of I/O with computation (our parser tries to do some computation while data is being 
read from the disk).
-- We have different versions of our parser. We use a heuristic to select the best one for the input 
file.
-- We also do some assumptions about the input format (e.g.: the global ids will always be in a format 
that is similar to the ones shown in the GISCUP website). If these assumptions do not hold true, an 
exception is thrown (when this is detected) and our main program selects another version of the algorithm 
(that makes fewer assumptions about the format, but should be slower…)

- We tested the algorithm with several input datasets we created (from small files with 10 features to 
files with several GBs and tens of millions of features).
- We also created test cases with several special cases (e.g.: disconnected graphs, graphs with multiple 
starting points/controllers, graphs without controllers, graphs without starting points, graphs with loops, 
graphs with multiple edges connecting the same pair of vertices, graphs with lines connecting several points, 
trees, datasets where the only starting point is also a controller, starting points being edges, graphs 
where a controller is an articulation point, graphs without edges where the only vertex is both a controller 
and a starting point, etc)

------------------------------------------------------------------------------------------------------
																	Description of the main files
-- main.cpp: selects which parser to use, calls it, calls the solver and writes the output file.
-- GiscupSolver.cpp: here we actually solve the problem, creating the BC-Tree and finding the vertices 
that are in a simple path between a controller and a starting point (the edges are found later based on 
the vertices).
-- GiscupGraph.h: defines the Graph we employ to represent the GISCUP problem.
-- GiscupParser.h: this is actually a base class that each parser must implement. In the base class we store 
the graph representing the problem and the information about the global id of each vertex/edge.
-- FastParallelGiscupParser.cpp: here we implement our fast custom parser.
-- RapidjsonParser.cpp: this is a parser based on the RapidJSON library (mentioned above). 
-- ImportantIdentifier.h: 
---- In this file we have the class ImportantIdentifier . It stores the “important” tokens found in the 
JSON file (e.g.: viaGlobalID, fromGlobalId, etc). These tokens are created by the parser and employed to create 
the graph. The tokens representing vertices are labeled with unique ids (0,1,...numVertices-1) and the ones 
representing edges are labeled similarly.
--- This file also defines the way the global IDs are stored (do we store their strings? Do we store pointers 
to the position of the ids in a big char array containing the entire JSON?  Do we encode the global ids using 
2 64-bit integers?)

