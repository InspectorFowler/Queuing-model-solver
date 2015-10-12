function [  ] = Graphu2powall( k,r,f,div)
t = cputime;
Graphu2(k,r,f,div,3)
Graphu2(k,r,f,div,2)
Graphu2(k,r,f,div,1)
Computation_time = cputime - t
end
