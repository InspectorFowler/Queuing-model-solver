function [ ] = CompSol3( k,r,f,un,pow )
Syms z x y w;
um = un/3;
uf = um/3;
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
    Am1 = subs(Am,z,1);
    Af = [g1 -(r*z) (z*un*x); -(f*z) g2 (z*um*y); 0 -(f*z) (z*uf*w)];
    Af = det(Af);
    Af1 = subs(Af,z,1);
    Q = factor(det(A))/(z-1);
    roots = double(solve(Q,z));
    real_roots = 0;
    imaginary_roots = 0;
    between_roots = 0;
    deg = size(roots,1);
    j=1;
    specroot = zeros(2,1);
    for i = 1:size(roots,1)
        if isreal(roots(i,1)) == 1
            real_roots = real_roots + 1;
            if roots(i,1)>0 && roots(i,1)<1
                between_roots = between_roots + 1;
                while j<3
                specroot(j,1)= roots(i,1);
                j=j+2;
                end
                j=j-1;
            end
        else isreal(roots(i,1)) == 0;
            imaginary_roots = imaginary_roots + 1;
        end
    end
    fprintf('\nDegree of polynomial Q(z): %d \n', deg)
    fprintf('Number of Real roots of Q(z): %d \n', real_roots)
    fprintf('Number of Imaginary roots of Q(z): %d \n', imaginary_roots)
    fprintf('Number of roots of Q(z) in (0,1): %d \n', between_roots)
    fprintf('The Roots being: ');roots  
    q = subs(Q,z,1);
    fprintf('Value of Q(1) = %d \n',q);
    %%----------------------------------------------------------------------
    %%Solution using An
    eq1n = subs(An,z,specroot(1,1));
    eq2n = subs(An,z,specroot(2,1));
    eq3n = ((q*gn)- An1);
    Psoln = solve(eq1n,eq2n,eq3n,x,y,w);
    Psoln = [double(Psoln.x); double(Psoln.y); double(Psoln.w)];
    fprintf('Value of Unknown probabilities: ');Psoln
    Gnet = (An + Am + Af)/Q;
    Gnet = subs(Gnet,x,Psoln(1,1));
    Gnet = subs(Gnet,y,Psoln(2,1));
    Gnet = subs(Gnet,w,Psoln(3,1));
    q_length = diff(Gnet,z);
    q_length = subs(q_length,z,1)/pow;
    fprintf('\nExpected number in the system : %f \n', q_length);
    %%----------------------------------------------------------------------
    %%Solution using Am
    eq1m = subs(Am,z,specroot(1,1));
    eq2m = subs(Am,z,specroot(2,1));
    eq3m = ((q*gm)- Am1);
    Psolm = solve(eq1m,eq2m,eq3m,x,y,w);
    Psolm = [double(Psolm.x); double(Psolm.y); double(Psolm.w)];
    fprintf('Value of Unknown probabilities: ');Psolm
    Gnet = (An + Am + Af)/Q;
    Gnet = subs(Gnet,x,Psoln(1,1));
    Gnet = subs(Gnet,y,Psoln(2,1));
    Gnet = subs(Gnet,w,Psoln(3,1));
    q_length = diff(Gnet,z);
    q_length = subs(q_length,z,1)/pow;
    fprintf('\nExpected number in the system : %f \n', q_length);
    %%----------------------------------------------------------------------
    %%Solution using Af
    eq1f = subs(Af,z,specroot(1,1));
    eq2f = subs(Af,z,specroot(2,1));
    eq3f = ((q*gf)- Af1);
    Psolf = solve(eq1f,eq2f,eq3f,x,y,w);
    Psolf = [double(Psolf.x); double(Psolf.y); double(Psolf.w)];
    fprintf('Value of Unknown probabilities: ');Psolf
    Gnet = (An + Am + Af)/Q;
    Gnet = subs(Gnet,x,Psoln(1,1));
    Gnet = subs(Gnet,y,Psoln(2,1));
    Gnet = subs(Gnet,w,Psoln(3,1));
    q_length = diff(Gnet,z);
    q_length = subs(q_length,z,1)/pow;
    fprintf('\nExpected number in the system : %f \n', q_length);
else
    fprintf('\nValue of the parameters do not satisfy the stability condition\n')
end
end

