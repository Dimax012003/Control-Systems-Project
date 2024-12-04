syms x1 x2;


V=0.2*x1^(2)-0.2*x1*x2+0.3*x2^(2);
V_dot=-0.2*x2^(2)-0.2*x1^(2)-0.2*x1^(2)*x2+0.5*x1*x2^(2)+0.3*x2^(3);

x=-0.5:0.01:0.5;
y=-0.5:0.01:0.5;

V_list=[];
x_sol(1,1)=0;
x_sol(1,2)=0;
k=1;
for i=[1:length(x)]
    for j=[1:length(y)]
        if double(subs(V_dot,{x1,x2},{double(x(i)),double(y(j))}))>0

            V_list=[V_list double(subs(V,{x1,x2},{double(x(i)),double(y(j))}))];
            x_sol(k,1)=x(i);
            x_sol(k,2)=y(j);
            k=k+1;
            break;
        end
    end
end

