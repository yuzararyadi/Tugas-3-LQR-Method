% simulate_quadcopter_LQR.m
% Simulasi LQR untuk subsystems attitude, position, altitude
% Berdasarkan Okyere et al. (2019). Referensi: file paper yang Anda upload. :contentReference[oaicite:6]{index=6}

clear; close all; clc;

%% 0) PARAMETERS (ganti sesuai data nyata jika ada)
m = 1.2;        % mass (kg) - typical small quadrotor
Ix = 0.014;     % inertia around x (kg*m^2)
Iy = 0.014;     % inertia around y
Iz = 0.028;     % inertia around z
g = 9.81;

% For linear position model we need effective 'a' = 1/(m*U1_hover). 
% We approximate U1_hover so that thrust = mg: U1h = m*g
U1h = m*g;
a = 1/(m*U1h);  % simplified linearization constant used in paper

%% 1) Attitude subsystem (states: phi, phi_dot, theta, theta_dot, psi, psi_dot)
A_att = [ 0 1 0 0 0 0;
          0 0 0 0 0 0;
          0 0 0 1 0 0;
          0 0 0 0 0 0;
          0 0 0 0 0 1;
          0 0 0 0 0 0];
% B as in paper maps U2,U3,U4 (we use 3 inputs for attitude)
B_att = [0 0 0;
         1/Ix 0 0;
         0 0 0;
         0 1/Iy 0;
         0 0 0;
         0 0 1/Iz];

%% 2) Position subsystem (x, x_dot, y, y_dot)
A_pos = [0 1 0 0;
         0 0 0 0;
         0 0 0 1;
         0 0 0 0];
% B_pos maps Ux and Uy (paper used 'a' factor)
B_pos = [0 0;
         a 0;
         0 0;
         0 a];

%% 3) Altitude subsystem (z, z_dot)
A_alt = [0 1; 0 0];
B_alt = [0; 1]; % in paper B = [0; 1] with Uz as input after normalization

%% 4) Initial Q,R using Bryson's rule (choose max acceptable values)
% Attitude: consider max angles 30 deg (~0.5236 rad) and angular rate 5 rad/s
max_phi = deg2rad(30); max_phidot = 5;
max_theta = deg2rad(30); max_thetadot = 5;
max_psi = deg2rad(180); max_psidot = 5;

Qatt = diag([1/max_phi^2, 1/max_phidot^2, 1/max_theta^2, 1/max_thetadot^2, 1/max_psi^2, 1/max_psidot^2]);
Ratt = 1*eye(3) * (1/(100^2)); % assume max motor input squared ~ placeholder

% Position: max pos error 5 m, max vel 3 m/s
Qpos = diag([1/5^2, 1/3^2, 1/5^2, 1/3^2]);
Rpos = 0.5*eye(2) * (1/(50^2));

% Altitude: max z error 5 m, zdot 3 m/s
Qalt = diag([1/5^2, 1/3^2]);
Ralt = 1 * (1/(100^2));

%% 5) Compute LQR gains
K_att = lqr(A_att, B_att, Qatt, Ratt);  % 3x6
K_pos = lqr(A_pos, B_pos, Qpos, Rpos);  % 2x4
K_alt = lqr(A_alt, B_alt, Qalt, Ralt);  % 1x2

disp('LQR gains:');
disp('K_att:'); disp(K_att);
disp('K_pos:'); disp(K_pos);
disp('K_alt:'); disp(K_alt);

%% 6) Build closed-loop A_cl for each subsystem (for linear simulation)
Acl_att = (A_att - B_att*K_att);
Acl_pos = (A_pos  - B_pos*K_pos);
Acl_alt = (A_alt  - B_alt*K_alt);

%% 7) Simulate step response (initial condition or reference)
T = 0:0.01:10; % 10 s simulation

% Example 1: attitude - command zero initial disturbance -> start with initial state
x0_att = [deg2rad(10); 0; deg2rad(5); 0; 0; 0]; % initial tilt disturbance
[y_att, t_att, x_att] = lsim(ss(Acl_att, B_att, eye(6), 0), zeros(length(T),3), T, x0_att);

% Compute control u(t) = -K * x(t) for attitude
u_att = -(K_att * x_att')'; % each row sample -> columns are inputs

% Example 2: position - command step reference via adding reference as input:
% For simplicity simulate initial position offset and let controller drive to zero
x0_pos = [1.5; 0; -0.8; 0]; % initial position errors
[y_pos, t_pos, x_pos] = lsim(ss(Acl_pos, B_pos, eye(4), 0), zeros(length(T),2), T, x0_pos);
u_pos = -(K_pos * x_pos')';

% Example 3: altitude
x0_alt = [0.8; 0]; % initial altitude error
[y_alt, t_alt, x_alt] = lsim(ss(Acl_alt, B_alt, eye(2), 0), zeros(length(T),1), T, x0_alt);
u_alt = -(K_alt * x_alt')';

%% 8) Plotting results (states + control signals)
figure('Name','Attitude states and control'); 
subplot(3,1,1); plot(T, rad2deg(x_att(:,[1 3 5]))); ylabel('angle (deg)');
legend('\phi','\theta','\psi'); grid on;
subplot(3,1,2); plot(T, x_att(:,[2 4 6])); ylabel('ang vel (rad/s)'); legend('\phi_dot','\theta_dot','\psi_dot'); grid on;
subplot(3,1,3); plot(T, u_att); ylabel('u (att inputs)'); xlabel('Time (s)'); legend('U2','U3','U4'); grid on;

figure('Name','Position states and control');
subplot(2,1,1); plot(T, x_pos(:,[1 3])); ylabel('pos (m)'); legend('x','y'); grid on;
subplot(2,1,2); plot(T, u_pos); ylabel('u (pos inputs)'); xlabel('Time (s)'); legend('Ux','Uy'); grid on;

figure('Name','Altitude state and control');
subplot(2,1,1); plot(T, x_alt(:,1)); ylabel('z (m)'); grid on;
subplot(2,1,2); plot(T, u_alt); ylabel('u_z'); xlabel('Time (s)'); grid on;

%% 9) Evaluate performance metrics (rise time, overshoot, settling time)
% We'll compute for attitude phi (x_att(:,1)) as example
phi = rad2deg(x_att(:,1));
t = T;
% target is 0 (we measured decay to 0). For rise time we consider time to go from 10%->90% of initial magnitude.
y0 = phi(1); y_final = phi(end);
[ymax, idx_max] = max(phi);
overshoot = (ymax - y0)/abs(y0) * 100; % percent relative to initial magnitude

% rise time (10%->90% of initial magnitude toward final)
threshold_low = 0.1*y0; threshold_high = 0.9*y0;
idx_low = find(abs(phi) <= abs(threshold_low),1); idx_high = find(abs(phi) <= abs(threshold_high),1);
if isempty(idx_low) || isempty(idx_high)
    rise_time = NaN;
else
    rise_time = t(idx_high) - t(idx_low);
end
% settling time: time after which |phi - y_final| < 2% of y0 for all t
tol = 0.02*abs(y0);
idx_settle = find(abs(phi - y_final) > tol);
if isempty(idx_settle)
    settling_time = 0;
else
    settling_time = t(max(idx_settle));
end

fprintf('Attitude phi metrics: overshoot=%.2f%%, rise_time=%.3f s, settling_time=%.3f s\n', ...
    overshoot, rise_time, settling_time);

%% END
