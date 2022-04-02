% Copyright (C) 2021 Xiaobo Zheng - All Rights Reserved
% developed by Xiaobo at BIT Jan 2021

%% Setup

clc; 
clear; 
close all;
% rng(0);
psi_f = 225/180*pi;

%% UAV problem definition
scale_para = 200;

% Initial state
x_0 = [5000.0/scale_para; 5000.0/scale_para; 5000.0/scale_para; 300.0; 0/180*pi; 180/180*pi];

% Target state
x_star = [0.0; 0.0; 0.0; 300; -60/180*pi; psi_f];

% Time horizon
t_f = 35.0;

% Time interval
dt = 0.2;

% Maximum magnitude of control
u_max = [30; 30];

% Initialize dynamics
fprintf("initializing inverted UAV dynamics...\n")
m = 150.0; g = 9.8;
dyn = UAV3D_Dynamics(m, g, scale_para);

% Initialize cost
fprintf("initializing quadratic cost function...\n")
%%%%% PN
% if psi_f == 225/180*pi
%     Q_f     = blkdiag( 5.0*scale_para*eye(3), 0.0, 5*10^4*eye(2)); % dt = 0.2
% elseif psi_f == 240/180*pi
%     Q_f     = blkdiag( 5.0*scale_para*eye(3), 0.0, 5*10^4*eye(2)); % dt = 0.2
% elseif psi_f == 260/180*pi
%     Q_f     = blkdiag( 12.0*scale_para*eye(3), 0.0, 5*10^4*eye(2)); % dt = 0.2
% end
%%%%% UC
if psi_f == 225/180*pi
    Q_f     = blkdiag( 2.0*scale_para*eye(3), 0.0, 5*10^4*eye(2)); % dt = 0.2
elseif psi_f == 240/180*pi
    Q_f     = blkdiag( 2.0*scale_para*eye(3), 0.0, 5*10^4*eye(2)); % dt = 0.2 % for a higher precise, 4.0*scale
elseif psi_f == 260/180*pi
    Q_f     = blkdiag( 2.0*scale_para*eye(3), 0.0, 5*10^4*eye(2)); % dt = 0.2
end
Q_s     = blkdiag( 0.0000001 *scale_para* eye(3), 0.0, 0.0000001 * eye(2));
R       = 0.1 * eye(numel(u_max));
cost	= QuadraticCost3d(Q_f, Q_s, R);

% % Initialize obstacle
% fprintf("initializing obstacle...\n")
constraints = [];
center = [1500/scale_para, 2000/scale_para]; r = 500/scale_para;
constraints = [constraints, CircleConstraints_3D(center, r)];
center = [2800/scale_para, 4800/scale_para]; r = 500/scale_para;
constraints = [constraints, CircleConstraints_3D(center, r)];

% Number of DDP iterations
num_iter = 300;

% DDP line search parameter
alpha = 0.5;

% Video framerate
fps = 30;

%% Execution of DDP

fprintf("executing DDP...");

tic;
sol = UAV3D_fft_cddp(x_0, x_star, t_f, dt, dyn, cost, constraints, u_max, num_iter, alpha);
toc;

%% Begin post-processing of solution

if sol.error == 1
    fprintf("Error\n");
end

% Extract the pendulum angle information from the solution
x = zeros(1, length(sol.x));
y = zeros(1, length(sol.x));
z = zeros(1, length(sol.x));
v = zeros(1, length(sol.x));
gama = zeros(1, length(sol.x));
psi = zeros(1, length(sol.x));
ay = zeros(1, length(sol.u));
az = zeros(1, length(sol.u));
for k = 1:length(sol.x)
    x(k) = sol.x{k}(1);
    y(k) = sol.x{k}(2);
    z(k) = sol.x{k}(3);
    v(k) = sol.x{k}(4);
    gama(k) = sol.x{k}(5);
    psi(k) = sol.x{k}(6);
end
for k = 1:length(sol.u)
    ay(k) = sol.u{k}(1);
    az(k) = sol.u{k}(2);
end


%% Plot trajectories

% Plot x history
figure;
pbaspect([5 3 1])
hold on;
plot(sol.t, x*scale_para, "k", "LineWidth", 2);
grid on;
xlabel("Time [s]", "Interpreter", "latex", "FontSize", 20);
ylabel("x  [m]", "Interpreter", "latex", "FontSize", 20);
ax = gca();
ax.FontSize = 16;
ax.TickLabelInterpreter = "latex";

% Plot y history
figure;
pbaspect([5 3 1])
hold on;
plot(sol.t, y*scale_para, "k", "LineWidth", 2);
grid on;
xlabel("Time [s]", "Interpreter", "latex", "FontSize", 20);
ylabel("y [m]", "Interpreter", "latex", "FontSize", 20);
ax = gca();
ax.FontSize = 16;
ax.TickLabelInterpreter = "latex";

% Plot z history
figure;
pbaspect([5 3 1])
hold on;
plot(sol.t, z*scale_para, "k", "LineWidth", 2);
grid on;
xlabel("Time [s]", "Interpreter", "latex", "FontSize", 20);
ylabel("z [m]", "Interpreter", "latex", "FontSize", 20);
ax = gca();
ax.FontSize = 16;
ax.TickLabelInterpreter = "latex";

% Plot v history
figure;
pbaspect([5 3 1])
hold on;
plot(sol.t, v, "k", "LineWidth", 2);
grid on;
xlabel("Time [s]", "Interpreter", "latex", "FontSize", 20);
ylabel("v [m]", "Interpreter", "latex", "FontSize", 20);
ax = gca();
ax.FontSize = 16;
ax.TickLabelInterpreter = "latex";

% Plot anglular gama history
figure;
pbaspect([5 3 1])
hold on;
plot(sol.t, rad2deg(gama), "k", "LineWidth", 2);
grid on;
xlabel("Time [s]", "Interpreter", "latex", "FontSize", 20);
ylabel("Anglular gama [rad]", "Interpreter", "latex", "FontSize", 20);
ax = gca();
ax.FontSize = 16;
ax.TickLabelInterpreter = "latex";

% Plot anglular psi history
figure;
pbaspect([5 3 1])
hold on;
plot(sol.t, rad2deg(psi), "k", "LineWidth", 2);
grid on;
xlabel("Time [s]", "Interpreter", "latex", "FontSize", 20);
ylabel("Anglular psi [rad]", "Interpreter", "latex", "FontSize", 20);
ax = gca();
ax.FontSize = 16;
ax.TickLabelInterpreter = "latex";

% Plot trajectory in state space
figure;
pbaspect([5 3 1])
hold on;
axis equal
% view(3)
plot3(x.*scale_para, y.*scale_para, z.*scale_para, "k", "LineWidth", 2);
plot3(x_0(1)*scale_para, x_0(2)*scale_para, x_0(3)*scale_para, "o", "MarkerFaceColor", "blue", ...
                          "MarkerEdgeColor", "blue");
plot3(x_star(1)*scale_para, x_star(2)*scale_para, x_star(3)*scale_para, "o", "MarkerFaceColor", ...
                                     "green", "MarkerEdgeColor", "green");
for i = 1:length(constraints)
    [hCylinder, hEndPlate1, hEndPlate2] = cylinder3([constraints(i).Center(1) constraints(i).Center(2) x_0(3)].*scale_para, ...
        [constraints(i).Center(1) constraints(i).Center(2) x_star(3)].*scale_para, constraints(i).r.*scale_para, 1000, 'r', 1, 0);
end
grid on;
xlabel("X [m]", "Interpreter", "latex", "FontSize", 20);
ylabel("Y [m]", "Interpreter", "latex", "FontSize", 20);
zlabel("Z [m]", "Interpreter", "latex", "FontSize", 20);
xlim([-500 6000]); ylim([0 6000]); zlim([0 6000]);
ax = gca();
ax.FontSize = 16;
ax.TickLabelInterpreter = "latex";

% Plot control sequence ay
figure;
pbaspect([5 3 1])
hold on;
plot(sol.t(1:end-1), ay, "b", "LineWidth", 2);
% plot(sol.t, ay, "b", "LineWidth", 2);
grid on;
xlabel("Time [s]", "Interpreter", "latex", "FontSize", 20);
ylabel("Control Input ay [N]", "Interpreter", "latex", "FontSize", 20);
ax = gca();
ax.FontSize = 16;
ax.TickLabelInterpreter = "latex";

% Plot control sequence az
figure;
pbaspect([5 3 1])
hold on;
plot(sol.t(1:end-1), az, "b", "LineWidth", 2);
% plot(sol.t, az, "b", "LineWidth", 2);
grid on;
xlabel("Time [s]", "Interpreter", "latex", "FontSize", 20);
ylabel("Control Input az [N]", "Interpreter", "latex", "FontSize", 20);
ax = gca();
ax.FontSize = 16;
ax.TickLabelInterpreter = "latex";

% Plot cost function vs iteration
figure;
pbaspect([5 3 1])
hold on;
plot(1:length(sol.J), sol.J, "r", "LineWidth", 2);
grid on;
xlabel("DDP Iteration [-]", "Interpreter", "latex", "FontSize", 20);
ylabel("Cost Function [-]", "Interpreter", "latex", "FontSize", 20);
ax = gca();
ax.FontSize = 16;
ax.TickLabelInterpreter = "latex";
ax.YScale = 'log';

% Plot control energy vs iteration
figure;
pbaspect([5 3 1])
hold on;
plot(1:length(sol.E), sol.E, "r", "LineWidth", 2);
grid on;
xlabel("DDP Iteration [-]", "Interpreter", "latex", "FontSize", 20);
ylabel("Control Energy Usage [$\rm{N}^{2}\rm{m}^{2}\rm{s}$]", "Interpreter", "latex", "FontSize", 20);
ax = gca();
ax.FontSize = 16;
ax.TickLabelInterpreter = "latex";
ax.YScale = 'log';

% Plot terminal time vs iteration
figure;
pbaspect([5 3 1])
hold on;
plot(1:length(sol.Tf), sol.Tf, "r", "LineWidth", 2);
grid on;
xlabel("DDP Iteration [-]", "Interpreter", "latex", "FontSize", 20);
ylabel("Terminal time [$\rm{N}^{2}\rm{m}^{2}\rm{s}$]", "Interpreter", "latex", "FontSize", 20);
ax = gca();
ax.FontSize = 16;
ax.TickLabelInterpreter = "latex";
ax.YScale = 'log';
