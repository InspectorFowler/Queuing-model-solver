function [ q_length ] = CompShort2cust( k,r,f,un,div,pow )
Syms z x y;
uf = un/div;
gn = r/(r+f);
gf = f/(r+f);
condition = (un*gn) + (uf*gf);
if condition > (pow*k)
    g1 = (k*z*(1-z^pow))+(f*z)+(un*(z-1));
    g2 = (k*z*(1-z^pow))+(r*z)+(uf*(z-1));
    A = [g1 -(r*z); -(f*z) g2];    
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
    n=10;
    q = subs(Q,z,1);
    eq1f = subs(Af,z,specroot);
    eq2f = ((q*(f/(r+f))) - Af1);
    Psolf = solve(eq1f,eq2f,x,y);
    Psolf = [double(Psolf.x); double(Psolf.y)];
    Probs = zeros(((n*pow)+1),2);    
    Probs(1,1) = Psolf(1,1);
    Probs(1,2) = Psolf(2,1);
    Probs(2,1) = ((((k+f)*Probs(1,1))-(r*Probs(1,2)))/un);
    Probs(2,2) = ((((k+r)*Probs(1,2))-(f*Probs(1,1)))/uf);    
    for i=1:(pow-1)
        Probs(i+2,1)= (((k+f+un)*Probs(i+1,1))-(r*Probs(i+1,2)))/un;
        Probs(i+2,2)= (((k+r+uf)*Probs(i+1,2))-(f*Probs(i+1,1)))/uf;
    end    
    for i= pow:((n*pow)-1)
        Probs(i+2,1)= ((((k+f+un)*Probs(i+1,1))-(r*Probs(i+1,2))-(k*Probs((i+1-pow),1)))/un);
        Probs(i+2,2)= ((((k+r+uf)*Probs(i+1,2))-(f*Probs(i+1,1))-(k*Probs((i+1-pow),2)))/uf);
    end     
    Custprobs = zeros((n+1),2);
    Custprobs(1,1)=Probs(1,1);
    Custprobs(1,2)=Probs(1,2);
    for i=0:(n-1)
        totaln = 0;
        totalf = 0;
        for j=((pow*i)+2):((pow*(i+1))+1)
            totaln = totaln+Probs(j,1);
            totalf = totaln+Probs(j,2);
        end
        Custprobs((i+2),1)=totaln;
        Custprobs((i+2),2)=totalf;
    end
    q_length = 0;    
    for i=0:n
        q_length= q_length+(i*((Custprobs(i+1,1))+(Custprobs(i+1,2))));
    end
  else
    fprintf('\nValue of the parameters do not satisfy the stability condition\n')
end
end
