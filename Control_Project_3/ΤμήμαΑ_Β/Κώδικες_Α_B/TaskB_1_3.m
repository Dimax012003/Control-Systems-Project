x10=[11];
x20=[1];
x30=[1];
for i = 1:length(x10)
    [t,x]=ode45(@odefun,[0 10],[x10(i);x20(i);x30]);
    figure();
    plot(t,x(:,1),t,x(:,2),t,x(:,3));
    title(['x_0= [',num2str(x10(i)),' ',num2str(x20(i)),']']);
    legend('x1','x2','θ_{hat}');
    xlabel('t');
    grid on;
    saveas(gcf,['graphs',num2str(i),'.png']);
end


function dx=odefun(t,x)
    dx=[-x(1)+x(2);-x(1)+(0.5-x(3))*x(2)^2;5*x(2)^3-2*x(2)^2*x(1)];
end