clc; clear; close all;

% Motor parameters (Phase III)
K_m = 0.1;
R_f = 2;
L_f = 0.5;
J   = 0.01;
b   = 0.001;

s = tf('s');

% Plant (speed transfer function)
G_omega = K_m / ((L_f*s + R_f)*(J*s + b));  % V_f -> omega

% Display open-loop info
disp('Open-loop poles:'), pole(G_omega)
disp(['Open-loop DC gain G_omega(0) = ', num2str(dcgain(G_omega))])

% PID candidates (manual)
candidates = {
    'Conservative', 1.0, 0.5, 0.05;
    'Moderate',     5.0, 1.0, 0.5;
    'Aggressive',   12.0, 2.0, 1.5
};

t = linspace(0,50,2000);

% Plot open-loop step
figure;
step(G_omega, t)
title('Open-loop step response (speed \omega(t))')
grid on
saveas(gcf, 'fig_openloop_step.png')

% Evaluate each controller
for k = 1:size(candidates,1)
    name = candidates{k,1};
    Kp = candidates{k,2};
    Ki = candidates{k,3};
    Kd = candidates{k,4};
    
    C = Kp + Ki/s + Kd*s;   % PID
    L = C * G_omega;
    T = feedback(L, 1);
    
    % Plot closed-loop step
    figure;
    step(T, t)
    title(sprintf('Closed-loop step: %s PID (Kp=%.2g, Ki=%.2g, Kd=%.2g)', name, Kp, Ki, Kd))
    grid on
    saveas(gcf, sprintf('fig_cl_step_%s.png', name));
    
    % Print poles and stepinfo
    p_cl = pole(T);
    disp(['--- ', name, ' ---']);
    disp('Closed-loop poles:'); disp(p_cl.');
    S = stepinfo(T);
    disp('Step info:'); disp(S);
    disp(['DC gain (closed-loop): ', num2str(dcgain(T))]);
end
