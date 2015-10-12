function [ q_length ] = CompShort3( k,r,f,un,pow )
um = un*0.5;
uf = un*0.25;
Syms z x y w;
gn = ((r/f)^2)/(1+(r/f)+((r/f)^2));
gm = (r/f)/(1+(r/f)+((r/f)^2));
gf = 1/(1+(r/f)+((r/f)^2));
condition = (un*gn) + (um*gm) + (uf*gf);
if condition > (pow*k)
    g1 = (k*z*(1-z^pow))+(f*z)+(un*(z-1));
    g2 = (k*z*(1-z^pow))+(r*z)+(f*z)+(um*(z-1));
    g3 = (k*z*(1-z^pow))+(r*z)+(uf*(z-1));
    A = [g1 -(r*z) 0; -(f*z) g2 -(r*z); 0 -(f*z) g3];
    An = [(z*un*x) -(r*z) 0; (z*um*y) g2 -(r*z); (z*uf*w) -(f*z) g3];
    An = det(An);
    An1 = subs(An,z,1);
    Am = [g1 (z*un*x) 0; -(f*z) (z*um*y) -(r*z); 0 (z*uf*w) g3];
    Am = det(Am);
    Af = [g1 -(r*z) (z*un*x); -(f*z) g2 (z*um*y); 0 -(f*z) (z*uf*w)];
    Af = det(Af);    
    Q = factor(det(A))/(z-1);
    roots = double(solve(Q,z));
    j=1;
    specroot = zeros(2,1);
    for i = 1:size(roots,1)
        if roots(i,1)>0 && roots(i,1)<1
            while j<3
                specroot(j,1)= roots(i,1);
                j=j+2;
            end
            j=j-1;
        end
    end
    q = subs(Q,z,1);
    eq1n = subs(An,z,specroot(1,1));
    eq2n = subs(An,z,specroot(2,1));
    eq3n = ((q*gn)- An1);
    Psoln = solve(eq1n,eq2n,eq3n,x,y,w);
    Psoln = [double(Psoln.x); double(Psoln.y); double(Psoln.w)];
    Gnet = (An + Am + Af)/Q;
    Gnet = subs(Gnet,x,Psoln(1,1));
    Gnet = subs(Gnet,y,Psoln(2,1));
    Gnet = subs(Gnet,w,Psoln(3,1));
    q_length = diff(Gnet,z);
    q_length = subs(q_length,z,1);    
else
    fprintf('\nValue of the parameters do not satisfy the stability condition\n')
end
end


