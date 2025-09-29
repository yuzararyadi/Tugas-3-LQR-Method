% run_lqr_compare.m
% Reproduce results from: "Optimizing LQR controllers: A comparative study"
% (paper provided). See file for reference. :contentReference[oaicite:5]{index=5}

clear; close all; clc;

%% 1) System parameters (from Table 1)
mp = 0.23;   % mass of pendulum [kg]
M  = 2.4;    % mass of cart [kg]
I  = 0.099;  % inertia [kg m^2]
l  = 0.36;   % length to COM [m]
g  = 9.81;
b  = 0.05;   % friction cart
rho = 0.005; % rotary damping

% Precomputed intermediate values and A,B from paper (they already substituted values)
A = [0      1       0       0;
    -0.0909 -2.6755  0    -0.0227;
     0      0       0       1;
     0.2273 31.2136  0     0.2652];

B = [0; 1.8182; 0; -4.5455];

C = [1 0 0 0; 0 0 1 0];
D = zeros(2,1);

sys = ss(A,B,C,D);

%% 2) LQR with various Q/R from the paper (ABC, Proposed, Neural)
% ABC result (paper): Q and R given
Q_ABC = diag([45.3407, 8.4808, 56.4716, 0.2844]);
R_ABC = 0.75;
K_ABC = lqr(A,B,Q_ABC,R_ABC);

% Proposed method (paper)
Q_prop = diag([1.6244, 1.2260, 29.6911, 19.3185]);
R_prop = 0.01;
K_prop = lqr(A,B,Q_prop,R_prop);

% Neural LQR (paper)
Q_nn = diag([2.3017, 8.5611, 0.0001, 0.0005]);
R_nn = 0.0218;
K_nn = lqr(A,B,Q_nn,R_nn);

% For reference, get poles
p_ABC  = eig(A - B*K_ABC);
p_prop = eig(A - B*K_prop);
p_nn   = eig(A - B*K_nn);

disp('Closed-loop poles (ABC):');  disp(p_ABC.');
disp('Closed-loop poles (Proposed):'); disp(p_prop.');
disp('Closed-loop poles (Neural):'); disp(p_nn.');

%% 3) Impulse response comparisons (pendulum angle theta and cart pos x)
T = 0:0.01:10; % simulation time
% state output y = [x; theta] per paper C matrix
sys_cl_ABC  = ss(A - B*K_ABC,  B*0, C, 0); % impulse: apply impulse to state? better to use state impulse
sys_cl_prop = ss(A - B*K_prop, B*0, C, 0);
sys_cl_nn   = ss(A - B*K_nn,  B*0, C, 0);

% Simulate impulse on state: simulate impulse response of closed-loop for an initial state
% We'll use initial condition: small initial angle (theta0)
x0 = [0; 0; 0.1; 0]; % small angular perturbation (0.1 rad)
[y_abc, t]  = initial(sys_cl_ABC, x0, T);
[y_prop, t] = initial(sys_cl_prop, x0, T);
[y_nn,  t]  = initial(sys_cl_nn,  x0, T);


% Plot theta (second output)
figure;
plot(t, y_abc(:,2),'LineWidth',1.2); hold on;
plot(t, y_prop(:,2),'LineWidth',1.2);
plot(t, y_nn(:,2),'LineWidth',1.2);
xlabel('Time (s)'); ylabel('\theta (rad)');
legend('ABC-LQR','Proposed','Neural-LQR');
title('Pendulum angle response (initial \theta = 0.1 rad)');
grid on;

% Plot cart position x (first output)
figure;
plot(t, y_abc(:,1),'LineWidth',1.2); hold on;
plot(t, y_prop(:,1),'LineWidth',1.2);
plot(t, y_nn(:,1),'LineWidth',1.2);
xlabel('Time (s)'); ylabel('x (m)');
legend('ABC-LQR','Proposed','Neural-LQR');
title('Cart position response (initial \theta = 0.1 rad)');
grid on;

%% 4) Compute basic metrics: overshoot, settling time (theta)
% small helper to compute overshoot and settling time
metrics = @(tt, y) struct(...
    'OS', (max(y)-y(1))/abs(y(1))*100, ...
    'settling_time', settling_time(tt,y,0.02)); % 2% band

m_abc = metrics(t, y_abc(:,2));
m_prop = metrics(t, y_prop(:,2));
m_nn = metrics(t, y_nn(:,2));

fprintf('\nMetrics (pendulum theta) - Overshoot (%%) and settling time (s):\n');
fprintf('ABC:  OS=%.2f, Ts=%.2f\n', m_abc.OS, m_abc.settling_time);
fprintf('Proposed: OS=%.2f, Ts=%.2f\n', m_prop.OS, m_prop.settling_time);
fprintf('Neural: OS=%.2f, Ts=%.2f\n', m_nn.OS, m_nn.settling_time);

%% 5) Hitung sinyal kontrol u(t) = -K * x(t)
% Note: x(t) bisa didapat dari fungsi 'initial' juga, output tambahan
[x_abc, t]  = initial(ss(A-B*K_ABC,B,C,D), x0, T);  % ambil output default
[~,~,x_abc_states] = initial(ss(A-B*K_ABC,B,C,D), x0, T);

[~,~,x_prop_states] = initial(ss(A-B*K_prop,B,C,D), x0, T);
[~,~,x_nn_states]   = initial(ss(A-B*K_nn,B,C,D),   x0, T);

u_abc  = - (K_ABC  * x_abc_states')';   % kontrol input ABC
u_prop = - (K_prop * x_prop_states')';  % kontrol input Proposed
u_nn   = - (K_nn   * x_nn_states')';    % kontrol input Neural

%% supporting function: simple settling time estimate
function Ts = settling_time(t,y,band)
    y_final = y(end);
    idx = find(abs(y - y_final) > band*abs(y_final), 1, 'last');
    if isempty(idx)
        Ts = 0;
    else
        Ts = t(min(length(t), idx+1));
    end
end

%% 5) Hitung sinyal kontrol u(t) = -K * x(t)
% Note: x(t) bisa didapat dari fungsi 'initial' juga, output tambahan
[x_abc, t]  = initial(ss(A-B*K_ABC,B,C,D), x0, T);  % ambil output default
[~,~,x_abc_states] = initial(ss(A-B*K_ABC,B,C,D), x0, T);

[~,~,x_prop_states] = initial(ss(A-B*K_prop,B,C,D), x0, T);
[~,~,x_nn_states]   = initial(ss(A-B*K_nn,B,C,D),   x0, T);

u_abc  = - (K_ABC  * x_abc_states')';   % kontrol input ABC
u_prop = - (K_prop * x_prop_states')';  % kontrol input Proposed
u_nn   = - (K_nn   * x_nn_states')';    % kontrol input Neural

% Plot control input
figure;
plot(t, u_abc,'LineWidth',1.2); hold on;
plot(t, u_prop,'LineWidth',1.2);
plot(t, u_nn,'LineWidth',1.2);
xlabel('Time (s)'); ylabel('Control Input u(t) [N]');
legend('ABC-LQR','Proposed','Neural-LQR');
title('Control Input vs Time');
grid on;
