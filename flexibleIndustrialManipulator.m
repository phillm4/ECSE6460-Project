%% Robust Control of a Flexible Manipulator
% 
% ECSE 6460 Multivariable Control - Final Project
% Kimberly Oakes & Mitchell Phillips
% Last Edited: April 25, 2017

clc, clear, close all;
s = tf('s');

%% Model Notes
%
% Inputs to the system:
%                       u  - motor torque
%                       wm - motor disturbance
%                       wp - tool disturbance
%
% Outputs of the system:
%                       qm - motor position
%                       Pa - tool acceleration (Pdd)
%

%% Nominal Parameter Values
%
% Note for linearized model, k1 = constant and Td = 0
% k1 is the nonlinear stiffness for the gear box
%

Jm=5e-3; Ja1=2e-3; Ja2=0.02; Ja3=0.02; % moment of intertia, [kg*m^2]

k1_low = 16.7; % k1 min, [Nm/rad]
k1_high = 100; % k1 max, [Nm/rad]
k1_nom = linspace(k1_low, k1_high, 5); % linear set of k1 values

k2 = 110; k3 = 80; % stiffness, [Nm/rad]
d1 = 0.08; d2 = 0.06; d3 = 0.08; % damping, [Nm*s/rad]
fm=6e-3; fa1=1e-3; fa2=1e-3; fa3=1e-3; % viscous friction, [Nm*s/rad]
n = 220; % gear ratio
l1 = 20e-3; l2 = 600e-3; l3 = 1530e-3; % link lengths, [mm]
Td = 0; %0.5e-3; % time delay, [s]

%% Linearized Model

J = diag([Jm, Ja1, Ja2, Ja3]);

D = [d1, -d1, 0, 0;
    -d1, (d1+d2), -d2, 0;
    0, -d2, (d2+d3), -d2;
    0, 0, -d3, d3];

F = diag([fm, fa1, fa2, fa3]);

% in paper, used k1_high for analysis
K_high = [k1_high, -k1_high, 0, 0;
    -k1_high, (k1_high+k2), -k2, 0;
    0, -k2, (k2+k3), -k3;
    0, 0, -k3, k3];

A = [zeros(4,4), eye(4);
    -inv(J)*K_high, -inv(J)*(D+F)];

% only concerned with inputs on first and last mass
% B = [zeros(4,4);
%     diag([1/Jm, 0, 0, 1/Ja3])];
B = [zeros(4,4);
    inv(J)];

E = (1/n) * [0, l1, l2, l3, zeros(1,4)];
C = [1, zeros(1,7);
    E];

%% State Space and Transfer Function Systems
%
% based off the maximum stiffness value for the motor, k1_high
%

arm_ss = ss(A,B,C,[]); % state space system with k1_high
arm_ss.StateName = {'motor pos.  (rad)';'joint 1 pos. (rad)';...
    'joint 2 pos. (rad)';'joint 3 pos. (rad)';...
    'motor vel.  (rad/s)';'joint 1 vel. (rad/s)';...
    'joint 2 vel. (rad/s)';'joint 3 vel. (rad/s)'};
arm_ss.InputName = {'u+wm','null 1','null 2','wp'};  
arm_ss.OutputName = {'qm';'P'};

%% Transfer function, motor torque to motor acceleration, nominal moderl
%
% recreating plant response. Different figure numbers depending on paper.
% Included additional plants with various k1 values.
% Some values dependent on previous created parameters.
%

figure(1)
grid on
hold on
for i = 1:length(k1_nom)
    
    % stiffness matrix for each k1 value
    K_nom = [k1_nom(i), -k1_nom(i), 0, 0;
        -k1_nom(i), (k1_nom(i)+k2), -k2, 0;
        0, -k2, (k2+k3), -k3;
        0, 0, -k3, k3];
    A_nom = [zeros(4,4), eye(4);
        -inv(J)*K_nom, -inv(J)*(D+F)];
    
    B_nom = B;
    
    % output concerened with on motor velocity
    C_nom = [0,0,0,0,1,0,0,0];
    
    % generate system model
    sys_ss_nom = ss(A_nom,B_nom,C_nom,[]);
    
    % multipliy by 's' to get acceleration for output
    sys_tf_nom = tf(sys_ss_nom(1,1))*s; 
    
    h = bodeplot(sys_tf_nom);
    setoptions(h,'FreqUnits','Hz','PhaseVisible','off');
    xlim([1, 100])
    
end

title('Transfer Function, Motor Torque to Motor Acceleration')
grid on
legend(['k_{1} = ',num2str(k1_nom(1))],...
    ['k_{1} = ',num2str(k1_nom(2))],...
    ['k_{1} = ',num2str(k1_nom(3))],...
    ['k_{1} = ',num2str(k1_nom(4))],...
    ['k_{1} = ',num2str(k1_nom(5))])
hold off

%% Open Loop Analysis

% Poles and Zeros

% Output Poles
[T_SISO, Po] = eig(arm_ss.A);
Po = diag(Po);
Yp = arm_ss.C * T_SISO;
fprintf(...
    'Poles \t\t\t Output pole vectors^T \t Output pole directions^T \n')
disp([Po Yp' ]);

% Input Pole
[Q, Pi] = eig(arm_ss.A');
Pi = diag(Pi);

Up1 = arm_ss.B' * Q;
Q1 = Q(:, [1:4,7:8,5:6]);
Up1 = Up1(:, [1:4,7:8,5:6]);
Pi1 = Pi([1:4,7:8,5:6],:);

Up2 = [arm_ss.b(:,1),arm_ss.B(:,4)]' * Q;
Q2 = Q(:, [1:4,7:8,5:6]);
Up2 = Up2(:, [1:4,7:8,5:6]);
Pi2 = Pi([1:4,7:8,5:6],:);

fprintf([...
    'Poles \t\t\t\t Input pole vectors^T \t\t\t\t'...
    'Input pole directions^T \n'])
disp([Pi1 Up1' ]);

% transmission zeros
z = tzero(arm_ss({'qm','P'},{'u+wm','wp'}));
fprintf('System zeros: \n')
disp(z)

z_qm = tzero(arm_ss({'qm'},{'u+wm'}));
fprintf('Transmission zeros for Motor Position Output, u and wm Input: \n')
disp(z_qm)

%% Frequency Response 

G_SISO = arm_ss({'qm'},{'u+wm'}); % SISO Plant

G_MIMO = [arm_ss({'qm'},{'u+wm'});...
     arm_ss({'P'},{'u+wm'})*s^2]; % MIMO Plant

figure(2)
bodemag(G_MIMO,'b');
title('Open Loop Response')

%% Mode Response

[V,D] = eig(arm_ss.A,'nobalance');
p = diag(D);
[p,ndx] = esort(p);
V = V(:,ndx);
systemp = arm_ss({'qm','P'},{'u+wm'});
figure(3)
for ii = 1:length(D)
    b = real(V(:,ii))/norm(real(V(:,ii)));
    set(systemp,'b',b,'d',[]);
    subplot(2,4,ii);
    impulse(systemp);
    title(['Mode ',num2str(D(ii,ii))])
end

%% Loop Shaping

% SISO
S_SISO=1/(1+G_SISO); % SISO Sensitivity
T_SISO=1-G_SISO; % SISO Complementary 

% Gain and Phase Margins for SISO
Marg_SISO = allmargin(G_SISO)

% Controller 1 Weights - qm
w1_1 = 1;
w2_1 = 100*tf(conv([1 10],[0.5227 3.266 1406]),conv([1 0],[1 5.808 2324]));

% plot for loop shaping of controller 1
figure(4)
bodemag(S_SISO, {1e-1,1e3})
hold on
bodemag(1/w2_1, {1e-1,1e3})
title('Sensitivity and Performance Weight - Controller 1')
legend('Sensitivity', 'Weighting Function')
hold off

% MIMO
S_MIMO=1/(1+G_MIMO(1,1)); % MIMO Sensitivity
T_MIMO=1-G_MIMO(1,1); % MIMO Complementary 

% Gain and Phase Margins for MIMO
Marg_MIMO = allmargin(G_MIMO(1,1)) %should match up with Marg_SISO

% controller 2 weights - qm, Pdd
w1_2 = 50;
w2_2 = [tf([1 3],[1 0]) 0;0 tf(0.2,conv([1 5],[1 5]))];

% plot for loop shaping of controller 2
figure(5)
bodemag(S_MIMO, {1e-1,1e3})
hold on
bodemag((1/(w2_1(1,1)*w2_2(1,1))), {1e-1,1e3})
title('Sensitivity and Performance Weight - Controller 2')
legend('Sensitivity', 'Weighting Function')
hold off

% Controller Gain - do for both controllers
[kinf1,cl1,gam1,info1]=ncfsyn(arm_ss({'qm'},{'u+wm'}),w1_1,w2_1); % ctrl 1
[kinf2,cl2,gam2,info2] = ncfsyn(G_MIMO,w1_2,w2_2); % ctrl 2

% plot of controller gains
figure(6)
bodemag(kinf1,{1e-1,1e3},'r')
hold on
bodemag(kinf2(1,1),{1e-1,1e3},'b')
grid on
title('Controller Gains')
legend({'$H_{\infty} (q_{m})$','H_{\infty}(q_{m}, \ddot{P})'},...
    'Interpreter','latex')
hold off

% Loop Gains for both Controllers
L1 = kinf1*G_SISO;
L2 = kinf2(1)*G_MIMO(1);

figure(7)
bodemag(L1,{1e-1,1e3},'r')
hold on
bodemag(L2,{1e-1,1e3},'b')
grid on
title('Loop Gains')
legend({'$H_{\infty} (q_{m})$','H_{infty}(q_{m}, ddot{P})'},...
    'Interpreter','latex')
hold off


%% Loop Shaping Results

S1 = 1/(1+L1);
S2 = 1/(1+L2);

figure(8)
bodemag(S1,{1e-1,1e3},'r')
hold on
bodemag(1/w2_1, {1e-1,1e3},'r--')
bodemag(S2,{1e-1,1e3},'b')
bodemag((1/(w2_1(1,1)*w2_2(1,1))), {1e-1,1e3},'b--')
bodemag(S_SISO, {1e-1,1e3},'k.')
grid on
title('Controller Sensitivity')
legend({'$S (q_{m})$','$w (q_{m})$',...
    'S (q_{m}, ddot{P})','w (q_{m}, ddot{P})','S'},...
    'Interpreter','latex','Location','southeast')
hold off

%% Closed loop Analysis
% 
% sys = (kinf1*G_SISO)/(1-(kinf1*G_SISO)); % closed loop wrt ref 
% Use positive feedback for Hinf controller
%

disp('Closed Loop Poles using Controller 1')
pole(cl1)
disp('Closed Loop Zeros using Controller 1')
zero(cl1(1,1))
figure(9);
pzmap(cl1)
title('Pole/Zero Map Controller 1')

disp('Closed Loop Poles using Controller 2')
pole(cl2)
disp('Closed Loop Zeros using Controller 2')
zero(cl2(1,1))
figure(10);
pzmap(cl2)
title('Pole/Zero Map Controller 2')

figure(11);
step(cl1(1,1))
hold on
step(cl2(1,1))
title('Step Response for Output q_m and Input (u+w)')
legend('Controller 1', 'Controller 2')
hold off

figure(12);
step(cl1(2,1))
hold on
step(cl2(2,1))
title('Step Response for Output P_{dd} and Input (u+w)')
legend('Controller 1', 'Controller 2')
hold off

%% Singular Value Plots
%
% Does not use arm_ss model
%
% Replication of Fig. 4 in, "Singular Joint Control of a Flexible
% Industrial Manipulator using H_inf Loop Shaping
%
% Contains two plots for singular values for the system from u to y and w
% to y for y = q_m (blue) and y = [q_m Pdd]' (red)
% 

% recreate System, some parameters based off previous values
% 'sig' subscript used to represent 'singular value'

A_sig = A;

% only concerned with inputs u and w which act on fifth and eighth state
B_sig = [zeros(4,4);diag([1/Jm, 0, 0, 1/Ja3])];

C_sig = C;

% create systems for singular values
sys_ss_sig = ss(A_sig,B_sig,C_sig,[]);
sys_tf_sig = tf(sys_ss_sig);

% plot of the singular values, u to y
figure(13)
sigmaplot(sys_tf_sig(1,1),{1e-1,1e3},'b')
hold on
sigmaplot([(sys_tf_sig(1,1)) ; (sys_tf_sig(2,1)*s^2)],{1e-1,1e3},'r')
ylim([-100, 100])
title('Singular Values from u to y')
xlabel('Frequency [rads/s]')
ylabel('Magnitude [dB]')
legend('y = q_m','y = [q_m Pdd]''')
grid on
hold off

% plot of singular values, w to y
figure(14)
sigmaplot([sys_tf_sig(1,1), sys_tf_sig(1,4)],{1e-1,1e3},'b')
hold on
sigmaplot([sys_tf_sig(1,1), sys_tf_sig(1,4) ;...
    sys_tf_sig(2,1)*s^2,sys_tf_sig(2,4)*s^2],{1e-1,1e3},'r')
ylim([-100, 100])
title('Singular Values from w to y')
xlabel('Frequency [rads/s]')
ylabel('Magnitude [dB]')
legend('y = q_m','y = [q_m Pdd]''')
grid on
hold off

%% References
%
% https://www.researchgate.net/publication/228715342_Robust_control_of_a_flexible_manipulator_arm_A_benchmark_problem
%