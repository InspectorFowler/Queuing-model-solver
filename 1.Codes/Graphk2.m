function [  ] = Graphk2( r,f,un,pow )
uf = un/2;
kmax = ((r*un)+(f*uf))/((r+f)*pow);
kuse = kmax-0.075;
dataset = zeros(20,2);
for i=1:20
    dataset(i,1) = kuse;
    dataset(i,2) = CompShort2(kuse,r,f,un,pow);
    kuse = kuse*0.85;
end
k_values = dataset(1:20,1);
k_values_smooth = smooth(k_values,5,'lowess');
Ex_values = dataset(1:20,2);
temp = rand;
if temp>=0 && temp<=0.2
    p1 = plot(k_values_smooth, Ex_values,'c');
else if temp>0.2 && temp<=0.4
        p1 = plot(k_values_smooth, Ex_values,'m');
    else if temp>0.4 && temp<=0.6
            p1 = plot(k_values_smooth, Ex_values,'g');
        else if temp>0.6 && temp<=0.8
                p1 = plot(k_values_smooth, Ex_values,'k');
            else if temp>0.8 && temp<=1
                    p1 = plot(k_values_smooth, Ex_values,'r');
            end
        end
    end
end
xlabel('Arrival rate');
ylabel('Expected number in the system');
grid on;
Hold on;
end

