x10=[9 11];
x20=[11 1];



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
    ylabel('x_2');
    xlabel('x_1');
    grid on;
    saveas(gcf,['phase',num2str(i),'.png']);
end


function dx=odefun(t,x)
    dx=[-x(1)+x(2);-4*x(1)-5*x(2)];
end