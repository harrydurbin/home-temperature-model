%% CE 295 - Energy Systems and Control
%   HW 2 : Model Learning for Smart Home Thermal Management
%   Oski Bear, SID 18681868
%   Prof. Moura

% ode_lsqff.m
% ODEs for the gradient parameter identification algorithm
% t         : time
% y         : states of adaptive law, consisting of
%   theta_h   : parameter estimate
%   P         : covariance matrix
% data      : input-output data used to feed algorithm
% bbeta     : foregetting factor

function y_dot = ode_lsqff(t,y,data,bbeta)

%% Parse states of adaptive law
theta_h = y(1:2);
P = reshape(y(3:end),[2,2]);

%% Parse Input Data
it = data(:,1);     % Time vector
iT = data(:,2);     % Indoor temp. vector
iT_A = data(:,3);   % Ambient temp. vector
iT_B = data(:,4);   % Boiler temp. vector

%% Interpolate data
T = interp1(it,iT,t);
T_A = interp1(it,iT_A,t);
T_B = interp1(it,iT_B,t);

%% Parametric model notation
% Samping time step
dt = 1;

% Compute Room temperature at NEXT time step
T_plus = interp1(it,iT,t+dt);

% Compute \dot{T} using forward difference in time
z = (T_plus - T) / dt;

% Assemble regressor vector, \phi
phi = [T_A - T; T_B - T];

%% Gradient Update Law
% Normalization signal
msq = 1 + phi'*phi;

% Normalized Estimation error: \epsilon = (z - \theta_h^T \phi)/m^2
epsilon = (z - theta_h' * phi)/msq;

% Update Law
theta_h_dot = P * epsilon * phi;

% Update Law for Matrix P
P_dot = bbeta*P - P*((phi*phi')/msq)*P;

%% Concatenate derivatives
y_dot = [theta_h_dot; reshape(P_dot,[4,1])];