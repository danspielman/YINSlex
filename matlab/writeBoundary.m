function writeBoundary(where, what, fn)
% function writeBoundary(where, what, fn)
%
% write the where,what boundary conditions to file fn
% for communicating with java
% will produce a zero-indexed output

nb = length(where);
where = where - 1;

h = fopen(fn,'w');
fprintf(h,'%u\n',nb);
for i = 1:nb,
    fprintf(h,'%u %f\n',where(i),what(i));
end
