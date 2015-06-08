function [v,command] = lexFromFiles(a,where,what)
% function [v,command] = lexFromFiles(a,where,what)
%
% This code computes the lex minimizer by saving the graph and
% conditions in files, and then running java CompLexMinimizer

graphfile = tempname;
condfile = tempname;
outfile = tempname;

writeGraph(a,graphfile);
writeBoundary(where,what,condfile);
command = ['java -cp out/YINSlex_jar/YINSlex.jar CompLexMinimizer ', ...
      graphfile,' ', condfile, ' ', outfile];

unix(command);

v = load(outfile);

delete(graphfile);
delete(condfile);
delete(outfile);
