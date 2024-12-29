x10=[0.8,-0.4];
x20=[0.8,1];

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
    title(['Φασικό πορτρετο για ','x_0= [',num2str(x10(i)),' ',num2str(x20(i)),']']);
    hold on;
    scatter(x(:,1),x(:,2));
    grid on;
    xlabel('x1');
    ylabel('x2');
    saveas(gcf,['phase.png',num2str(i),'.png']);
end


function dx=odefun(t,x)
    P=[1 -0.5;-0.5 1.5];
    b=[0;1];
    dx=[-x(1)+x(2);-x(1)+x(1)*x(2)+0.5*x(2)^2-2*sign([x(1) x(2)]*P*b)*(x(1)^2+x(2)^2)];
end