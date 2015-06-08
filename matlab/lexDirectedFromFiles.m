function [v,command] = lexDirectedFromFiles(a,where,what)
% function [v,command] = lexDirectedFromFiles(a,where,what)
%
% This code computes the lex minimizer by saving the graph and
% conditions in files, and then running java CompLexDirected

graphfile = tempname;
condfile = tempname;
outfile = tempname;

writeDirectedGraph(a,graphfile);
writeBoundary(where,what,condfile);
command = ['java -cp out/YINSlex_jar/YINSlex.jar CompLexDirected ', ...
      graphfile,' ', condfile, ' ', outfile];

unix(command);

v = load(outfile);

delete(graphfile);
delete(condfile);
delete(outfile);
