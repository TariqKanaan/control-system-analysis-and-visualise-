%% ================================================================
% DC MOTOR ANALYSIS: THEORETICAL vs SIMULATION RESULTS
% Author: DR. Tariq (adapted by ChatGPT)
% ---------------------------------------------------------------
% Computes motor transfer functions, poles, DC gains, step response,
% impulse response, bode plots, and compares theoretical vs numerical.
% ================================================================

clc; clear; close all;

%% ---------------------------------------------------------------
% 1. MOTOR PARAMETERS
% ---------------------------------------------------------------
R_f = 2;          % Field resistance (ohm)
L_f = 0.5;        % Field inductance (H)
K_m = 0.1;        % Torque constant (Nm/A)
J   = 0.01;       % Rotor inertia (kg*m^2)
b   = 0.001;      % Viscous damping (N*m*s)

%% ---------------------------------------------------------------
% 2. THEORETICAL ANALYSIS
% ---------------------------------------------------------------
% Speed Transfer Function:
%   G_omega(s) = K_m / [(R_f + L_f s)(J s + b)]

den_omega_theoretical = [L_f*J, (L_f*b + R_f*J), (R_f*b)];
num_omega_theoretical = K_m;

% Poles (theoretical)
p_theoretical = roots(den_omega_theoretical);

% DC Gain (theoretical)
DC_gain_theoretical = K_m / (R_f * b);

%% ---------------------------------------------------------------
% 3. SIMULATION TRANSFER FUNCTIONS (MATLAB)
% ---------------------------------------------------------------
s = tf('s');

G_omega = K_m / ((R_f + L_f*s)*(J*s + b));
G_theta = G_omega / s;    % Integrator for position

% Poles (simulation)
p_sim = pole(G_omega);

% DC gain (simulation)
DC_gain_sim = dcgain(G_omega);

%% ---------------------------------------------------------------
% 4. COMPARISON TABLE (THEORETICAL vs SIMULATION)
% ---------------------------------------------------------------
fprintf('\n=== COMPARISON RESULTS ===\n');
fprintf('Theoretical Poles:\n');
disp(p_theoretical)

fprintf('Simulation Poles:\n');
disp(p_sim)

fprintf('Theoretical DC Gain: %f\n', DC_gain_theoretical);
fprintf('Simulation DC Gain : %f\n\n', DC_gain_sim);

%% ---------------------------------------------------------------
% 5. STEP RESPONSE (Speed and Position)
% ---------------------------------------------------------------
figure;
subplot(2,1,1)
step(G_omega);
title('Step Response - Speed \omega(t)');
grid on;

subplot(2,1,2)
step(G_theta);
title('Step Response - Position \theta(t)');
grid on;

%% ---------------------------------------------------------------
% 6. IMPULSE RESPONSE
% ---------------------------------------------------------------
figure;
subplot(2,1,1)
impulse(G_omega);
title('Impulse Response - Speed \omega(t)');
grid on;

subplot(2,1,2)
impulse(G_theta);
title('Impulse Response - Position \theta(t)');
grid on;

%% ---------------------------------------------------------------
% 7. BODE ANALYSIS
% ---------------------------------------------------------------
figure;
bode(G_omega);
grid on;
title('Bode Plot - Speed Transfer Function G_\omega(s)');

%% ---------------------------------------------------------------
% 8. ROOT LOCUS (OPTIONAL)
% ---------------------------------------------------------------
figure;
rlocus(G_omega);
title('Root Locus - Speed Transfer Function');

%% ---------------------------------------------------------------
% 9. PRINT SUMMARY
% ---------------------------------------------------------------
fprintf("Analysis complete.\n");
fprintf("Plots generated: step, impulse, bode, root locus.\n");
