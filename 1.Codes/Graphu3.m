function [  ] = Graphu3( k,r,f,pow )
gn = ((r/f)^2)/(1+(r/f)+((r/f)^2));
gm = (r/f)/(1+(r/f)+((r/f)^2));
gf = 1/(1+(r/f)+((r/f)^2));
umin = (pow*k)/(gn+(0.5*gm)+(0.25*gf));
uuse = umin+0.55;
dataset = zeros(20,2);
for i=1:20
    dataset(i,1) = uuse;
    dataset(i,2) = CompShort3(k,r,f,uuse,pow);
    uuse = uuse*1.1;
end
u_values = dataset(1:20,1);
u_values_smooth = smooth(u_values,5,'lowess');
Ex_values = dataset(1:20,2);
temp = rand;
if temp>=0 && temp<=0.2
    p1 = plot(u_values_smooth, Ex_values,'c');
else if temp>0.2 && temp<=0.4
        p1 = plot(u_values_smooth, Ex_values,'y');
    else if temp>0.4 && temp<=0.6
            p1 = plot(u_values_smooth, Ex_values,'g');
        else if temp>0.6 && temp<=0.8
                p1 = plot(u_values_smooth, Ex_values,'k');
            else if temp>0.8 && temp<=1
                    p1 = plot(u_values_smooth, Ex_values,'r');
            end
        end
    end
end
xlabel('Service Rate');
ylabel('Expected number in the system');
grid on;
Hold on;
end