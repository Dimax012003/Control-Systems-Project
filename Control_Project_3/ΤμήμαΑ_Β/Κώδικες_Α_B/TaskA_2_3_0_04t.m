r=zeros(1,6);
r_dot=0.04.*ones(1,6);


y=[-2,1,0,2,2.5,1.1];
y_dot=[0,0,0.5,2,-1,2];

e=r-y;
e_dot=r_dot-y_dot;


for i = 1:length(e_dot)
    [t,x]=ode45(@odefun1,[0 30],[e(i);e_dot(i)]);
    figure();
    rIN=0.04.*t;
    title(['x_0= [',num2str(e(i)),' ',num2str(e_dot(i)),']']);
    y_output=rIN-x(:,1);
    plot(t,x(:,1),t,x(:,2),t,y_output,t,rIN);
    legend('x1','x2','y','r');
    xlabel('t');
    grid on;
    saveas(gcf,['graphs',num2str(i),'.png']);

    figure();
    title(['x_0= [',num2str(e(i)),' ',num2str(e_dot(i)),']']);
    plot(x(:,1),x(:,2));
    hold on;
    scatter(x(:,1),x(:,2));
    grid on;
    xlabel('x1');
    ylabel('x2');
    saveas(gcf,['phase.png',num2str(i),'.png']);
end


function dx=odefun(t,x)
K=5;
T=0.2;
dx=[x(2);-(1/T)*x(2)-(K/T)*x(1)+1.2/T];
end

function dx=odefun1(t,x)
    K=5;
    T=0.2;
    a=0.05;
    if abs(x(1))<=0.2
        dx=[x(2);-(1/T)*x(2)-(K*a/T)*x(1)+0.04/T];
    elseif abs(x(1))>0.2
        dx=[x(2);-(1/T)*x(2)-(K/T)*x(1)+0.04/T];
    end

end

