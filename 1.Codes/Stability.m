function [  ] = Stability( k,r,f,un,div,pow )
Syms z x y;
uf = un/div;
gn = r/(r+f);
gf = f/(r+f);
condition1 = (pow*k)/(1/((gn/un) + (gf/uf)))
condition2 = (pow*k)/((un*gn) + (uf*gf))
end