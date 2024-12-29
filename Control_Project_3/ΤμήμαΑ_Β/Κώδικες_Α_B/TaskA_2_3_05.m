r=0.5.*ones(1,6);
r_dot=zeros(1,6);

y=[-2,1,0,2,2.5,1.1];
y_dot=[0,0,0.5,2,-1,2];

e=r-y;
e_dot=r_dot-y_dot;


for i = 1:length(e_dot)
    [t,x]=ode45(@odefun,[0 5],[e(i);e_dot(i)]);
    figure();
    rIN=0.5*ones(length(x(:,1)),1);

    y_output=rIN-x(:,1);
    plot(t,x(:,1),t,x(:,2),t,y_output,t,rIN);
    title(['x_0= [',num2str(e(i)),' ',num2str(e_dot(i)),']']);
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
dx=[x(2);-(1/T)*x(2)-(K/T)*x(1)];
end