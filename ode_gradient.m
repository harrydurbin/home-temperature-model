%% CE 295 - Energy Systems and Control
%   HW 2 : Model Learning for Smart Home Thermal Management
%   Harrison Durbin, SID 26951511
%   Prof. Moura

% ode_gradient.m
% ODEs for the gradient parameter identification algorithm
% t         : time
% theta_h   : parameter estimate
% data      : input-output data used to feed algorithm
% Gam       : Update law gain

function theta_h_dot = ode_gradient(t,theta_h,data,Gam)


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
% z = \dot{T} = (T(t+dt) - T(t))/dt
z = (T_plus-T)/dt ;

% Assemble regressor vector, \phi
phi = [T_A-T,T_B-T]' ;

%% Gradient Update Law
% Normalization signal
msq = 1+(phi'*phi) ;

% Estimation error: \epsilon = (z - \theta_h^T \phi /) msq
epsilon = (z-theta_h'*phi)/msq ;

% Update Law
theta_h_dot = Gam*phi*epsilon ; %Gam*epsilon*phi ;