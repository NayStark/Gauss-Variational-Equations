%% Gauss' Variational Equations for an Orbit

clc;
clear;
close all;

a_r = 8E-7;
a_t = 2E-7;
a_n = -3E-6;

mu=3.986E5;
a0= 15000;
e0= 0.5;
i0= deg2rad(28.5);
capOmega0= deg2rad(142);
omega0= deg2rad(74);
theta0= deg2rad(0);

ICs=[a0; e0; i0; capOmega0; omega0; theta0];

% For a given time span, from time 0s to time tf s.
tf = 2.592E6;

tspan1 = [0 tf];
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t_1, y_1] = ode45(@(t, y) GVE(t, y, mu, a_r, a_t, a_n), tspan1, ICs, options);

fin=y_1(end,:);
fin = [fin(1),fin(2),rad2deg(fin(3:6))];

fprintf('\t----PART A----\nFinal Orbital Elements after 30 days:\n');
fprintf('a = %.6f km\n', fin(1));
fprintf('e = %.6f\n', fin(2));
fprintf('i = %.6f deg\n', fin(3));
fprintf('Omega = %.6f deg\n', fin(4));
fprintf('omega = %.6f deg\n', fin(5));
fprintf('theta = %.6f deg\n', fin(6));

% For one orbital period
T_period = 2*pi*sqrt(a0^3/mu);

tspan2 = [0 T_period];

[t_2, y_2] = ode45(@(t, y) GVE(t, y, mu, a_r, a_t, a_n), tspan2, ICs, options);


y_2 = y_2;
y_2(:, 3:6) = rad2deg(y_2(:, 3:6));


delta_a = y_2(:,1) - ICs(1);
delta_e = y_2(:,2) - ICs(2);
delta_i = y_2(:,3) - ICs(3);
delta_Omega = y_2(:,4) - ICs(4);
delta_omega = y_2(:,5) - ICs(5);
delta_theta = y_2(:,6) - ICs(6);

figure;
subplot(3,2,1);
plot(t_2, delta_a, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('\Deltaa (km)');
title('Change in Semi-major Axis vs Time');
grid on;

subplot(3,2,2);
plot(t_2, delta_e, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('\Deltae');
title('Change in Eccentricity vs Time');
grid on;

subplot(3,2,3);
plot(t_2, delta_i, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('\Deltai (deg)');
title('Change in Inclination vs Time');
grid on;

subplot(3,2,4);
plot(t_2, delta_Omega, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('\Delta\Omega (deg)');
title('Change in RAAN vs Time');
grid on;

subplot(3,2,5);
plot(t_2, delta_omega, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('\Delta\omega (deg)');
title('Change in Argument of Perigee vs Time');
grid on;

subplot(3,2,6);
plot(t_2, delta_theta, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('\Delta\theta (deg)');
title('Change in True Anomaly vs Time');
grid on;


function sol = GVE(t, y, mu, a_r, a_t, a_n)
    
    a= y(1);
    e= y(2);
    i= y(3);
    capOmega = y(4);
    omega = y(5);
    theta = y(6);
    
    p=a*(1-e^2);
    h=sqrt(mu*p);
    r = p/(1 + e*cos(theta));

    da = ((2*a^2)/h) * (e*sin(theta)*a_r + p/r * a_t);
    de = 1/h * (p*sin(theta)*a_r + ((p+r)*cos(theta) + r*e)*a_t);
    di = (r*cos(theta+omega))/h * a_n;
    dcapOmega = (r*sin(theta+omega))/(h*sin(i)) * a_n;
    domega = 1/(h*e) * (-p*cos(theta)*a_r - (p+r)*sin(theta)*a_t);
    dtheta = h/r^2 + 1/(e*h) * (p*cos(theta)*a_r - (p+r)*sin(theta)*a_t);

    sol = [da; de; di; dcapOmega; domega; dtheta];

end