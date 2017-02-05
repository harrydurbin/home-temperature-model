%% CE 295 - Energy Systems and Control
%   HW 2 : Model Learning for Smart Home Thermal Management
%   Harrison Durbin, SID 26951511
%   Prof. Moura

% DURBIN_HARRISON_HW2.m

clear; close all;
fs = 15;    % Font Size for plots

%% Load Data
data = csvread('HW2_Data.csv'); 
t = data(:,1);      % t   : time vector [min]
thr = data(:,2);    % thr : time vector [hours]
T = data(:,3);      % T   : Home interior temperature [deg C]
T_A = data(:,4);    % T_A : Ambient outdoor temperature [deg C]
T_B = data(:,5);    % T_B : Boiler temperature [deg C]

%% Problem 2(b) - Persistence of Excitation

%%%% OPTION with 2-D parameter vector
phi = [T_A-T,T_B-T]'; % signals for 2D parametric model
t_end = t(end);
PE_mat = zeros(2);

phi_sq = zeros(2,2,length(t));
for k = 1:length(t)
    phi_sq(:,:,k) = phi(:,k) * phi(:,k)';
end
PE_mat(1,1) = 1/t_end * trapz(t, phi_sq(1,1,:));
PE_mat(2,1) = 1/t_end * trapz(t, phi_sq(2,1,:));
PE_mat(1,2) = 1/t_end * trapz(t, phi_sq(1,2,:)); 
PE_mat(2,2) = 1/t_end * trapz(t, phi_sq(2,2,:)); 

PE_lam_min = min(eig(PE_mat)); %  MINIMUM EIGENVALUE OF PE_mat
fprintf(1,'PE Level for 2D Version : %1.4f\n',PE_lam_min);

%%%% OPTION with 3-D parameter vector
phi = [T_A,T_B,T]'; % signals for 3D parametric model
t_end = t(end);
PE_mat = zeros(3);

phi_sq = zeros(3,3,length(t));
for k = 1:length(t)
    phi_sq(:,:,k) = phi(:,k) * phi(:,k)';
end
PE_mat(1,1) = 1/t_end * trapz(t, phi_sq(1,1,:));
PE_mat(2,1) = 1/t_end * trapz(t, phi_sq(2,1,:));
PE_mat(3,1) = 1/t_end * trapz(t, phi_sq(3,1,:));
PE_mat(1,2) = 1/t_end * trapz(t, phi_sq(1,2,:));
PE_mat(2,2) = 1/t_end * trapz(t, phi_sq(2,2,:));
PE_mat(3,2) = 1/t_end * trapz(t, phi_sq(3,2,:));
PE_mat(1,3) = 1/t_end * trapz(t, phi_sq(1,3,:));
PE_mat(2,3) = 1/t_end * trapz(t, phi_sq(2,3,:));
PE_mat(3,3) = 1/t_end * trapz(t, phi_sq(3,3,:));

PE_lam_min = min(eig(PE_mat)); % MINIMUM EIGENVALUE OF PE_mat
fprintf(1,'PE Level for 3D Version : %1.4f\n',PE_lam_min);

%% Problem 4(a)
% Plot T(t), T_A(t), T_B(t) versus 24 hours of time
figure(1)
plot(thr,T,'r-',thr,T_A,'b--',thr,T_B,'m:','LineWidth',2);
xlim([0 24])
ylim([0 35])
xlabel('Time [hr]')
ylabel('Temperature [deg C]','FontSize',fs), set(gca,'FontSize',fs);
legend('T, Home Interior','T_A, Environment','T_B, Boiler',fs);
grid on
grid minor

%% Problem 4(b)
% Assemble Data
data = [t, T, T_A, T_B];

% Initial conditions
theta_hat0 = [0.1,0.1]';  

% Update Law Gain
Gam_a = 1e-5*eye(2) ; % increase from super small
Gam_b = 1e-4*eye(2) ; 
Gam_c = 1e-3*eye(2) ; 
Gam_d = 1e-2*eye(2) ;   
Gam_e = 1e-1*eye(2) ;  

% Integrate ODEs
[~,y_a] = ode23s(@(t,y_a) ode_gradient(t,y_a,data,Gam_a), t, theta_hat0); % output y, which is theta_hat
[~,y_b] = ode23s(@(t,y_b) ode_gradient(t,y_b,data,Gam_b), t, theta_hat0);
[~,y_c] = ode23s(@(t,y_c) ode_gradient(t,y_c,data,Gam_c), t, theta_hat0);
[~,y_d] = ode23s(@(t,y_d) ode_gradient(t,y_d,data,Gam_d), t, theta_hat0);
[~,y_e] = ode23s(@(t,y_e) ode_gradient(t,y_e,data,Gam_e), t, theta_hat0);

% Parse output
theta_hat_a = y_a;
theta_hat_b = y_b;
theta_hat_c = y_c;
theta_hat_d = y_d;
theta_hat_e = y_e;

%% Problem 4(c)

theta_hat_1_a = theta_hat_a(:,1);
theta_hat_2_a = theta_hat_a(:,2);

theta_hat_1_b = theta_hat_b(:,1);
theta_hat_2_b = theta_hat_b(:,2);

theta_hat_1_c = theta_hat_c(:,1);
theta_hat_2_c = theta_hat_c(:,2);

theta_hat_1_d = theta_hat_d(:,1);
theta_hat_2_d = theta_hat_d(:,2);

theta_hat_1_e = theta_hat_e(:,1);
theta_hat_2_e = theta_hat_e(:,2);

% Plot parameter estimates
figure(2); clf;

plot(thr,theta_hat_1_e,'r-',thr,theta_hat_2_e,'b--','LineWidth',3);
hold on
plot(thr,theta_hat_1_a,'-',thr,theta_hat_2_a,'--','LineWidth',0.5);
plot(thr,theta_hat_1_b,'-',thr,theta_hat_2_b,'--','LineWidth',0.5);
plot(thr,theta_hat_1_c,'-',thr,theta_hat_2_c,'--','LineWidth',0.5);
plot(thr,theta_hat_1_d,'-',thr,theta_hat_2_d,'--','LineWidth',0.5);
xlim([0 24])
ylim([.04 0.14])
xlabel('Time [hr]')
ylabel('Parameter Estimates (Gradient Method)','FontSize',fs), set(gca,'FontSize',fs);
legend('\theta_h_1 (\Gamma=1e-1 x I)','\theta_h_2 (\Gamma=1e-1 x I)',fs);
grid on
grid minor

% Output exit estimates
theta_hat = theta_hat_e ;
fprintf(1,'theta_hat1 (Gradient Method) : %1.4f\n',theta_hat(end,1));
fprintf(1,'theta_hat2 (Gradient Method) : %1.4f\n',theta_hat(end,2));

% True model parameters
R1 = 2 ;                 % Thermal resistance w/ outisde [deg C/kW]
R2 = 0.75 ;              % Thermal resistance w/ boiler [deg C/kW]
C = 10  ;                 % Thermal capacitance [kWh/deg C]
theta_1_true = [1/(C*R1)] ;
theta_2_true = [1/(C*R2)] ;

% Output true parameter values
fprintf(1,'theta_true1 : %1.4f\n',theta_1_true);
fprintf(1,'theta_true2 : %1.4f\n',theta_2_true);

%% Problem 5(a)
% System matrices for identified model
Ahat =  [-theta_hat(end,1)-theta_hat(end,2)]; % must be nxn matrix -> 1x1
Bhat =  [theta_hat(end,1),theta_hat(end,2)]; % 1x2 matrix

% Output states only (dummy variables, not used later)
C_dummy = 1; % 1x1 matrix
D_dummy = [[0],[0]]; % 2x1 matrix

% State space model
sys_hat = ss(Ahat,Bhat,C_dummy,D_dummy);

%% Problem 5(b)
% Load Validation Data
valdata = csvread('HW2_valData.csv');
t=valdata(:,1);             % t   : time vector [min]
thr=valdata(:,2);   % thr : time vector [hours]
Tval=valdata(:,3); % Tval   : Home interior temperature [deg C]
T_Aval=valdata(:,4); % T_Aval : Ambient outdoor temperature [deg C]
T_Bval=valdata(:,5); % T_Bval : Boiler temperature [deg C]

% Input vector from validation data set
U_hat =  [[T_Aval],[T_Bval]] ;

% Initial condition for indoor temp. [deg C]
That0 = 22  ;

% Simulate
[~,~,That] = lsim(sys_hat, U_hat, t, That0);

% Plot predicted and actual indoor temperature from validation data set
figure(3); clf;
plot(thr,Tval,'r-',thr,That,'b--','Linewidth',2)
xlim([0 24])
ylim([15 25])
xlabel('Time [hr]')
ylabel('Temperature [deg C]','FontSize',fs), set(gca,'FontSize',fs);
legend('Tval, From Validation Data','That, Estimated w/ Gradient Method',fs);
grid on
grid minor


%% Problem 6 - Least Squares with Forgetting Factor

% Forgetting factor
beta_ = .01;

% Initial conditions
theta_hat0 = [0.1,0.1]';  
P_0 = 1e-1*eye(2); % assume initially is equal to Gamma
P_0_reshaped = reshape(P_0, 4, 1);
theta_hat0__P_0 = [theta_hat0 ; P_0_reshaped];

% Integrate ODEs
[~,y_conc] = ode23s(@(t,y_conc) ode_lsq(t,y_conc,data,beta_), t, theta_hat0__P_0);

% Parse output
theta_hat_ = y_conc(:,1:2);

theta_hat_1_ = theta_hat_(:,1);
theta_hat_2_ = theta_hat_(:,2);

% Output exit estimates
fprintf(1,'theta_hat1 (LSQ Method): %1.4f\n',theta_hat_(end,1));
fprintf(1,'theta_hat2 (LSQ Method): %1.4f\n',theta_hat_(end,2));

% Plot parameter estimates
figure(4); clf;
plot(thr,theta_hat_1_,'r-',thr,theta_hat_2_,'b--','Linewidth',3);
xlim([0 24])
ylim([.04 0.14])
xlabel('Time [hr]')
ylabel('Estimated Parameters (LSQ w/FF Method)','FontSize',fs), set(gca,'FontSize',fs);
legend('\theta_h_1 (Beta=0.99)','\theta_h_2 (Beta=0.99)',fs);
grid on
grid minor

% System matrices for identified model
Ahat =  [-theta_hat_(end,1)-theta_hat_(end,2)]; % must be nxn matrix -> 1x1 
Bhat =  [theta_hat_(end,1),theta_hat_(end,2)]; % 1x2 matrix


% Output states only (dummy variables, not used later)
C_dummy = 1; % 1x1 matrix
D_dummy = [[0],[0]]; % 2x1 matrix

% State space model
sys_hat = ss(Ahat,Bhat,C_dummy,D_dummy);
 

% Load Validation Data
valdata = csvread('HW2_valData.csv');
t=valdata(:,1);             % t   : time vector [min]
thr=valdata(:,2);   % thr : time vector [hours]
Tval=valdata(:,3); % Tval   : Home interior temperature [deg C]
T_Aval=valdata(:,4); % T_Aval : Ambient outdoor temperature [deg C]
T_Bval=valdata(:,5); % T_Bval : Boiler temperature [deg C]

% Input vector from validation data set
U_hat =  [[T_Aval],[T_Bval]] ;

% Initial condition for indoor temp. [deg C]
That0 = 22  ;

% Simulate
[~,~,That] = lsim(sys_hat, U_hat, t, That0);

% Plot predicted and actual indoor temperature from validation data set
figure(5); clf;
plot(thr,Tval,'r-',thr,That,'b--','Linewidth',2)
xlim([0 24])
ylim([15 25])
xlabel('Time [hr]')
ylabel('Temperature [deg C]','FontSize',fs), set(gca,'FontSize',fs);
legend('Tval, From Validation Data','That, Estimated w/LSQ w/FF Method',fs);
grid on
grid minor
