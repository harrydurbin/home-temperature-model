%% CE 295 - Energy Systems and Control
%   HW 2 : Model Learning for Smart Home Thermal Management
%   Harrison Durbin, SID 26951511
%   Prof. Moura

% ode_lsq.m
% ODEs for the gradient parameter identification algorithm
% t         : time
% y_conc    : the concatenate of theta_h and P, used for running the two ODEs
% simultaneously
% theta_h   : parameter estimate
% P         : LSQ factor
% data      : input-output data used to feed algorithm
% beta      : forgetting factor


function theta_h_dot__P_dot = ode_lsq(t,y_conc,data,beta)


%% Parse Input Data
it = data(:,1);     % Time vector
it_end=it(end);     % End of time vector
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

%% Least Squares Algorithm (LSQ) w/ forgetting factor

theta_h = y_conc(1:2,1);
P=reshape(y_conc(3:6),2,2);

% Normalization signal
msq = 1+(phi'*phi) ;

% Estimation error: \epsilon = z - \theta_h^T \phi
epsilon = (z-(theta_h'*phi))/msq  ;

% Update Law
theta_h_dot = P*epsilon*phi;   
P_dot = beta*P-P*(phi'*phi)/msq*P;

P_dot_reshaped = reshape(P_dot, 4, 1);
theta_h_dot__P_dot=[theta_h_dot;P_dot_reshaped];


