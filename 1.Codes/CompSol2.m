function [ ] = CompSol2( k,r,f,un,pow )
Syms z x y;
uf = un/2;
gn = r/(r+f);
gf = f/(r+f);
condition = (un*gn) + (uf*gf);
kmin = ((r*un)+(f*uf))/((r+f)*pow);
if condition > (pow*k)
    g1 = (k*z*(1-z^pow))+(f*z)+(un*(z-1));
    g2 = (k*z*(1-z^pow))+(r*z)+(uf*(z-1));
    A = [g1 -(r*z); -(f*z) g2];
    An = [(z*un*x) -(r*z); (z*uf*y) g2];
    An = det(An);
    An1 = subs(An,z,1);
    Af = [g1 (z*un*x); -(f*z) (z*uf*y)];
    Af = det(Af);
    Af1 = subs(Af,z,1);
    Q = factor(det(A))/(z-1);
    roots = double(solve(Q,z));
    real_roots = 0;
    imaginary_roots = 0;
    between_roots = 0;
    deg = size(roots,1);
    for i = 1:size(roots,1)
        if isreal(roots(i,1)) == 1
            real_roots = real_roots + 1;
            if roots(i,1)>0 && roots(i,1)<1
                between_roots = between_roots + 1;
                specroot = roots(i,1);
            end
        else isreal(roots(i,1)) == 0;
            imaginary_roots = imaginary_roots + 1;
        end
    end
    fprintf('\nDegree of polynomial Q(z): %d \n', deg)
    fprintf('Number of Real roots of Q(z): %d \n', real_roots)
    fprintf('Number of Imaginary roots of Q(z): %d \n', imaginary_roots)
    fprintf('Number of roots of Q(z) in (0,1): %d \n', between_roots)
    fprintf('Value of Q(1) =  \n');(un*r)+(uf*f)-(pow*k*(r+f))
    q = (un*r)+(uf*f)-(pow*k*(r+f));
    fprintf('The Roots being: ');roots  
    %%------------------------------------------------------------------
    %%solution using An
    eq1n = subs(An,z,specroot);
    eq2n = ((q*(r/(r+f))) - An1);
    Psoln = solve(eq1n,eq2n,x,y);
    Psoln = [double(Psoln.x); double(Psoln.y)];
    fprintf('Value of Unknown probabilities: ');Psoln
    Gnet = (An + Af)/Q;
    Gnet = subs(Gnet,x,Psoln(1,1));
    Gnet = subs(Gnet,y,Psoln(2,1));
    q_length = diff(Gnet,z);
    q_length = subs(q_length,z,1)/pow;
    fprintf('\nExpected number in the system : %f \n', q_length);
    %%------------------------------------------------------------------
    %%solution using Af
    eq1f = subs(Af,z,specroot);
    eq2f = ((q*(f/(r+f))) - Af1);
    Psolf = solve(eq1f,eq2f,x,y);
    Psolf = [double(Psolf.x); double(Psolf.y)];
    fprintf('\nValue of Unknown probabilities: ');Psolf
    Gnet = (An + Af)/Q;
    Gnet = subs(Gnet,x,Psolf(1,1));
    Gnet = subs(Gnet,y,Psolf(2,1));
    q_length = diff(Gnet,z);
    q_length = subs(q_length,z,1)/pow;
    fprintf('\nExpected number in the system : %f \n', q_length);
  else
    fprintf('\nValue of the parameters do not satisfy the stability condition\n')
end
end
