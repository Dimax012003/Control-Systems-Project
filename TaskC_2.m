opts = odeset('Refine',1,'RelTol', 0.0000001, 'AbsTol', 0.00000001, 'MaxStep',0.1);
tspan = [0 20];
y0 = [pi/3 pi/3 0 0];
[t, x] = ode15s(@odefun, tspan, y0, opts);
xd1=pi/4+pi/6*sin(0.2*pi*t);
xd2=-pi/3+pi/3*cos(0.2*pi*t);

xd1_dot=pi^2/30*cos(0.2*pi*t);
xd2_dot=-pi^2/15*sin(0.2*pi*t);

figure();
plot(x(:,1)-xd1,x(:,3)-xd1_dot);
hold on;
scatter(x(:,1)-xd1,x(:,3)-xd1_dot);
hold off;
title('e_1  , de_1/dt');
xlabel('e_1');
ylabel('de_1/dt');
grid on;


figure();
plot(x(:,2)-xd2,x(:,4)-xd2_dot);
hold on;
scatter(x(:,2)-xd2,x(:,4)-xd2_dot);
hold off;
title('e_2  , de_2/dt');
xlabel('e_2');
ylabel('de_2/dt');
grid on;

figure();
plot(t, x(:,1), 'b-', 'LineWidth', 2); 
hold on;
plot(t, x(:,2), 'r-', 'LineWidth', 2);  
plot(t, xd1, 'g--', 'LineWidth', 2);   
plot(t, xd2, 'c--', 'LineWidth', 2);  
hold off;
xlabel('t(s)', 'FontSize', 12); 
ylabel('Aπόκριση θέσης ως προς τον χρόνο', 'FontSize', 12);
legend({'q_1', 'q_2', 'q_1d', 'q_2d'}, 'FontSize', 10, 'Location', 'best');
grid on;
set(gca, 'LineWidth', 1, 'FontSize', 12);

figure();
plot(t, x(:,3), 'b-', 'LineWidth', 2); 
hold on;
plot(t, x(:,4), 'r-', 'LineWidth', 2);  
plot(t, xd1_dot, 'g--', 'LineWidth', 2);   
plot(t, xd2_dot, 'c--', 'LineWidth', 2); 
hold off;
xlabel('t(s)', 'FontSize', 12); 
ylabel('Aπόκριση ταχύτητας ως προς τον χρόνο', 'FontSize', 12);
legend({'dq_1/dt', 'dq_2/dt','dq_1d/dt', 'dq_2d/dt'}, 'FontSize', 10, 'Location', 'best');
grid on;
set(gca, 'LineWidth', 1, 'FontSize', 12);  



figure();
plot(t, x(:,1)-xd1,'Color','blue','LineWidth',1);
hold on;
plot(t, x(:,2)-xd2,'Color','red','LineWidth',1);

xlabel('t(s)', 'FontSize', 12); 
ylabel('Σφάλμα θέσης ως προς τον χρόνο', 'FontSize', 12);
hold off;
legend('e_1', 'e_2');
grid on;

set(gca, 'LineWidth', 1, 'FontSize', 12);


function dxdt = odefun(t, x)

    m1 = 6;
    m2 = 4;
    ml = 0.5;
    lc1 = 0.2;
    lc2 = 0.1;
    l1 = 0.5;
    l2 = 0.4;
    I1 = 0.43;
    I2 = 0.05;
    g = 9.81;

    % Υπολογισμός του H
    h11 = m1*lc1^2 + m2*(lc2^2 + l1^2 + 2*l1*lc2*cos(x(2))) + ...
          ml*(l2^2 + l1^2 + 2*l1*l2*cos(x(2))) + I1 + I2;
    h12 = m2*lc2*(lc2 + l1*cos(x(2))) + ml*l2*(l2 + l1*cos(x(2))) + I2;
    h22 = lc2^2*m2 + l2^2*ml + I2;

    H = [h11, h12; h12, h22];

    % Υπολογισμός του C
    k = -l1*(m2*lc2 + ml*l2);

    c11 = k*sin(x(2))*x(4);
    c22 = 0;
    c12 = k*sin(x(2))*(x(3) + x(4));
    c21 = -k*sin(x(2))*x(3);

    C = [c11, c12; c21, c22];

    % Υπολογισμός του G
    g1 = (m2*lc2 + ml*l2)*g*cos(x(1) + x(2)) + (m2*l1 + ml*l1 + m1*lc1)*g*cos(x(1));
    g2 = (m2*lc2 + ml*l2)*g*cos(x(1) + x(2));

    g_w = [g1; g2];

    lc1_hat=0.25;
    lc2_hat=0.175;
    I1_hat=0.26;
    I2_hat=0.08;
    ml_hat=1;
    

    g1_hat = (m2*lc2_hat + ml_hat*l2)*g*cos(x(1) + x(2)) + (m2*l1 + ml_hat*l1 + m1*lc1_hat)*g*cos(x(1));
    g2_hat = (m2*lc2_hat + ml_hat*l2)*g*cos(x(1) + x(2));

    g_hat=[g1_hat;g2_hat];

    k_hat = -l1*(m2*lc2_hat + ml_hat*l2);
    c11_hat = k_hat*sin(x(2))*x(4);
    c22_hat = 0;
    c12_hat = k_hat*sin(x(2))*(x(3) + x(4));
    c21_hat = -k_hat*sin(x(2))*x(3);

    C_hat=[c11_hat,c12_hat;c21_hat,c22_hat];


    h11_hat = m1*lc1_hat^2 + m2*(lc2_hat^2 + l1^2 + 2*l1*lc2_hat*cos(x(2))) + ...
          ml_hat*(l2^2 + l1^2 + 2*l1*l2*cos(x(2))) + I1_hat + I2_hat;
    h12_hat = m2*lc2_hat*(lc2_hat + l1*cos(x(2))) + ml_hat*l2*(l2 + l1*cos(x(2))) + I2_hat;
    h22_hat = lc2_hat^2*m2 + l2^2*ml_hat + I2_hat;

    H_hat=[h11_hat,h12_hat;h12_hat,h22_hat];

    f=-inv(H)*(C*[x(3);x(4)]+g_w);
    G=inv(H);

    f_hat=-inv(H_hat)*(C_hat*[x(3);x(4)]+g_hat);
    G_hat=inv(H_hat);

    xd=[pi/4+pi/6*sin(0.2*pi*t);-pi/3+pi/3*cos(0.2*pi*t)];
    xd_dot=[pi^2/30*cos(0.2*pi*t);-pi^2/15*sin(0.2*pi*t)];
    xd_Ddot=[-pi^3/150*sin(0.2*pi*t);-pi^3/75*cos(0.2*pi*t)];
    d=29.35*sqrt(x(3)^2+x(4)^2)+24.2289+3.3168*2.935+(abs(x(3))+abs(x(4)))*norm(xd_dot+10*[x(1)-xd(1),x(2)-xd(2)]) + 0.55*(abs(x(3))+abs(x(4)))*sqrt((x(3))^2+(x(4))^2);

    
    ueq=-inv(G_hat)*(10*[x(3)-xd_dot(1);x(4)-xd_dot(2)]+f_hat-xd_Ddot);
   
    if(x(3)-xd_dot(1)+10*(x(1)-xd(1)))<=1e-5 &&  (x(3)-xd_dot(1)+10*(x(1)-xd(1)))>=-1e-5
        u(1)=ueq(1)-d*100000*(x(3)-xd_dot(1)+10*(x(1)-xd(1)));
    elseif (x(3)-xd_dot(1)+10*(x(1)-xd(1)))>1e-5
        u(1)=ueq(1)-d;
    elseif (x(3)-xd_dot(1)+10*(x(1)-xd(1)))<-1e-5
        u(1)=ueq(1)+d;
    end

    if(x(4)-xd_dot(2)+10*(x(2)-xd(2)))<=1e-5 && (x(4)-xd_dot(2)+10*(x(2)-xd(2)))>=-1e-5
        u(2)=ueq(2)-d*100000*(x(4)-xd_dot(2)+10*(x(2)-xd(2)));
    elseif (x(4)-xd_dot(2)+10*(x(2)-xd(2)))>1e-5
        u(2)=ueq(2)-d;
    elseif (x(4)-xd_dot(2)+10*(x(2)-xd(2)))<-1e-5
        u(2)=ueq(2)+d;
    end

    dxdt = zeros(4, 1);
    dxdt(1:2) = x(3:4);
    dxdt(3:4) = -H \ (C * x(3:4) + g_w)+G*u';
    fprintf('t: %.2f, d: %.2f, u: [%.2f %.2f]\n', t, d, u(1), u(2));
end