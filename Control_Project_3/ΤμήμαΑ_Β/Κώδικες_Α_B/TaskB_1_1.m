x10=[0.1 0.4 0.8 -0.4];
x20=[0.1 0.4 0.8 1];

x1=[-1:0.01:1];
x2=[-1:0.01:1];

syms x1 x2
eqn = 0.2*x1^2 - 0.3*x1*x2 + 0.2*x2^2 - 0.01057;

for i = 1:length(x10)
    [t,x]=ode45(@odefun,[0 10],[x10(i);x20(i)]);
    figure();
    plot(t,x(:,1),t,x(:,2));
    title(['x_0= [',num2str(x10(i)),' ',num2str(x20(i)),']']);
    legend('x1','x2');
    xlabel('t');
    grid on;
    saveas(gcf,['graphs',num2str(i),'.png']);


    figure();
    plot(x(:,1),x(:,2));
    title(['Φασικό πορτραιτο x_1,x_2 με ','x_0= [',num2str(x10(i)),' ',num2str(x20(i)),']']);
    hold on;
    scatter(x(:,1),x(:,2));
    hold on;
    fimplicit(eqn, [-1 1 -1 1],LineWidth=1); 
    ylabel('x_2');
    xlabel('x_1');
    legend('','','Πεδίο έλξης');
    grid on;
    saveas(gcf,['phase',num2str(i),'.png']);
end


function dx=odefun(t,x)
    dx=[-x(1)+x(2);-x(1)+x(1)*x(2)+0.5*x(2)^2];
end