function v = compLex(a,where,what)
% function v = compLex(a,where,what)
%
% the weights in a give the lengths of the edges: longer -> weaker connection
% a is the adjacency matrix of the input graph

import lexvolts.*;

[ai,aj,av] = find(a);
ijv = EdgeListGraph(ai-1,aj-1,av);
rlv = RandLexVolt(ijv,randi([0,1e10]));
rlv.setTerminals(where-1,what);
rlv.computeVoltages5;
v = rlv.voltages;
