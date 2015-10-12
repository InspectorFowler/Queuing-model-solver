function [ q_wait ] = CompShort2wait_alt( k,r,f,un,div,pow )
uf = un/div;
gn = r/(r+f);
gf = f/(r+f);
condition = (un*gn) + (uf*gf);
if condition > (pow*k)
    q_length = ((k*pow)/((gn*un+gf*uf)-(k*pow)))*((pow+1)/2)
    q_wait = ((k*pow)/((gn*un+gf*uf)-(k*pow)))*((pow+1)/2)*((gn/un)+(gf/uf));
    
else
    fprintf('\nValue of the parameters do not satisfy the stability condition\n')
end
end