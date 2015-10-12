function [ q_length ] = CompShort2( k,r,f,un,div,pow )
Syms z x y;
uf = un/div;
gn = r/(r+f);
gf = f/(r+f);
condition = (un*gn) + (uf*gf);
if condition > (pow*k)
    g1 = (k*z*(1-z^pow))+(f*z)+(un*(z-1));
    g2 = (k*z*(1-z^pow))+(r*z)+(uf*(z-1));
    A = [g1 -(r*z); -(f*z) g2];
    An = [(z*un*x) -(r*z); (z*uf*y) g2];
    An = det(An);    
    Af = [g1 (z*un*x); -(f*z) (z*uf*y)];
    Af = det(Af);
    Af1 = subs(Af,z,1);
    Q = factor(det(A))/(z-1);
    roots = double(solve(Q,z));
    for i = 1:size(roots,1)
        if roots(i,1)>0 && roots(i,1)<1
           specroot = roots(i,1);
        elseif roots(i,1)==0
               specroot = 0; 
        end
    end
    q = subs(Q,z,1);
    eq1f = subs(Af,z,specroot);
    eq2f = ((q*(f/(r+f))) - Af1);
    Psolf = solve(eq1f,eq2f,x,y);
    Psolf = [double(Psolf.x); double(Psolf.y)];
    Gnet = (An + Af)/Q;
    Gnet = subs(Gnet,x,Psolf(1,1));
    Gnet = subs(Gnet,y,Psolf(2,1));
    q_length = diff(Gnet,z);
    q_length = subs(q_length,z,1); 
  else
    fprintf('\nValue of the parameters do not satisfy the stability condition\n')
end
end
