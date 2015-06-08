function v = compLexDirected(a,where,what)
% function v = compLexDirected(a,where,what)
%
% the weights in a give the lengths of the edges: longer -> weaker connection
% a is the adjacency matrix

import lexvoltsdirected.*;

[ai,aj,av] = find(a);
rlv = lexvoltsdirected.RandLexVoltDirected(ai-1,aj-1,av,randi(1e5));
rlv.setTerminals(where-1,what);
rlv.computeVoltages();

n = length(a);
v = zeros(n,1);
vert= rlv.vertices;

 for i = 1:n
 v(i) = vert(i).voltage;
 end
 
