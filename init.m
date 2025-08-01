%% ===============================
%% Simulation Settings
%% ===============================

T_stop = 4 * 3600;    % [s] Total simulation time: 4 hours
dt      = 0.05;          % [s] Discretization time step (for numerical integration)

%% ===============================
%% Tank Geometry
%% ===============================

A1 = 1.5;             % [m^2] Cross-sectional area of Tank 1
A2 = 1.5;             % [m^2] Cross-sectional area of Tank 2
A3 = 1.5;             % [m^2] Cross-sectional area of Tank 3

%% ===============================
%% Flow Dynamics
%% ===============================

k1 = 0.8;             % [√(m^3/s)] Flow coefficient from Tank 1 → Tank 2
k2 = 0.6;             % [√(m^3/s)] Flow coefficient from Tank 2 → Tank 3
k_drain = 0.03;       % [1/s] Recycle coefficient from Tank 3 → Tank 1 (triggered by PID)

%% ===============================
%% Initial Conditions
%% ===============================

h1_0 = 0.4;           % [m] Initial water level in Tank 1
h2_0 = 0.3;           % [m] Initial water level in Tank 2
h3_0 = 0.2;           % [m] Initial water level in Tank 3

%% ===============================
%% Control Setpoints
%% ===============================

T1.setpoint = 0.5;    % [m] Target water level in Tank 1 (controlled by PID)
T3.setpoint = 0.4;    % [m] Threshold level in Tank 3 to trigger recycling (controlled)
