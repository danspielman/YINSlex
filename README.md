YINSlex
=======

This contains the code used in the experiments from the paper
  "Algorithms for Lipschitz Learning on Graphs" by Rasmus Kyng, Anup Rao, Sushant Sachdeva and Daniel Spielman,
  COLT 2015.  
  The paper is also available at http://arxiv.org/abs/1505.00290
  
This repository presently contains code for computing lexicographic minimizers in undirected and directed graphs.
The code is written in java, and is meant to be called either from Matlab or from the command line.

A jar file has been included, so it should not be necessary to compile the code.

To use the code from Matlab, type "init" in the main directory.
The two main routines in Matlab are compLex and compLexDirected.
compLex takes as input the adjacency matrix of a graph, followed by two vectors where and what giving the fixed voltages.  Where should be a vector of integers between 1 and the number of vertices, indicating which vertices have fixed values.  What indicates the values to which they are fixed.
Here is an example with a 6-by-6 grid graph.

```
>> a = grid2(6);
>> where = [1 6]; what = [0 1];
>> v = compLex(a,where,what);
>> reshape(v,6,6)

ans =

         0    0.1429    0.2449    0.3178    0.3698    0.4070
    0.2000    0.2857    0.3469    0.3907    0.4219    0.4442
    0.4000    0.4286    0.4490    0.4636    0.4740    0.4814
    0.6000    0.5714    0.5510    0.5364    0.5260    0.5186
    0.8000    0.7143    0.6531    0.6093    0.5781    0.5558
    1.0000    0.8571    0.7551    0.6822    0.6302    0.5930
```

It is also possible to call the code directly from the command line. 
In Unix-like environments, one can even do this from Matlab using the function
lexFromFiles.
It has the advantage of allowing the java VM more memory than Matlab usually provides.

To invoke the code directly from the command line, one must provide it with a graph,
a description of the boundary conditions, and the name of the file to which it should write the output.

An example of a graph is in data/grid3.txt.
The first line contains two integers: the number of vertices and the number of edges.
Every successive line contains one edge given as two integers indicating the vertices followed
by a float giving the length of the edge.
Note that vertices are indexed from 0, and longer lengths indicate looser connections.

An example of a file giving the boundary conditions is in data/boundary.txt.
The first line consists of one integer indicating how many vertices are on the boundary.
Each succesive line consists of an integer indicating the number of a vertex followed by its real number voltage.

You may run the code like:

```
java -cp out/YINSlex_jar/YINSlex.jar CompLexMinimizer data/grid3.txt data/boundary.txt voltages.txt
```

The output file, voltages.txt, contains one float on each line giving the voltage
at the corresponding vertex.

There are also versions of this code for directed graphs.
The matlab routines are compLexDirected and lexDirectedFromFiles.
The command line interface is CompLexDirected.

The most important thing to know about CompLexDirected is that in the listing
of edges in its input graph, the destination vertex comes before the source vertex.
Here is an example of the output of the code on a graph with 5 vertices and the voltages fixed at the endpoints.

```
$ cat data/path.txt
5 4
1 0 1.0
2 1 1.0
3 2 1.0
4 3 1.0

$ cat data/pathBdry.txt 
2
0 1
4 0
$ java -cp out/YINSlex_jar/YINSlex.jar CompLexDirected data/path.txt data/pathBdry.txt voltages.txt
$ cat voltages.txt 
1.0
0.75
0.5
0.25
0.0
$ cat data/pathBdry2.txt 
2
0 0
4 1
$ java -cp out/YINSlex_jar/YINSlex.jar CompLexDirected data/path.txt data/pathBdry2.txt voltages.txt
$ cat voltages.txt 
0.0
1.0
1.0
1.0
1.0
```

