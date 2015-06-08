function writeDirectedGraph(a, fn)
% function writeDirectedGraph(a, fn)
%
% write the undirected graph in a to the text file fn
% for communicating with java
% will produce a zero-indexed output

nv = length(a);
[ai,aj,av] = find(a);
ne = length(ai);

ai = ai-1;
aj = aj-1;

h = fopen(fn,'w');
fprintf(h,'%u %u\n',nv,ne);
for i = 1:ne,
    fprintf(h,'%u %u %f\n',ai(i),aj(i),av(i));
end
