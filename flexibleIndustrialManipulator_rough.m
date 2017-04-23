%% Robust Control of a Flexible Manipulator
% 
% ECSE 6460 Multivariable Control - Final Project
% Kimberly Oakes & Mitchell Phillips
% Last Edited: April 17, 2017

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
%                       Pa - tool acceleration
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
k1_nom = reshape(k1_nom, [1, 1, length(k1_nom)]);
z = ones(1,1,length(k1_nom)); % extend in z direction to create many models

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

K = [k1_nom, -k1_nom, 0*z, 0*z;
    -k1_nom, (k1_nom+k2*z), -k2*z, 0*z;
    0*z, -k2*z, (k2*z+k3*z), -k3*z;
    0*z, 0*z, -k3*z, k3*z];

A = [zeros(4,4), eye(4);
    -inv(J)*K, -inv(J)*(D+F)];

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
arm_ss.StateName = {'motor pos. (rad)';'joint 1 pos. (rad)';...
    'joint 2 pos. (rad)';'joint 3 pos. (rad)';...
    'motor vel. (rad/s)';'joint 1 vel. (rad/s)';...
    'joint 2 vel. (rad/s)';'joint 3 vel. (rad/s)'};
%sys_tf= tf(sys_ss); % transfer function system with k1_high

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

%% Singular Value Plots
%
% Replication of Fig. 4 in, "Singular Joint Control of a Flexible
% Industrial Manipulator using H_inf Loop Shaping
%
% Contains two plots for singular values for the system from u to y and w
% to y for y = q_m (blue) and y = [q_m Pdd]' (red)
% 

% recreate System, some parameters based off previous values
% 'sig' subscript used to represent 'singular value'

K_sig = K;

A_sig = A;

% only concerned with inputs u and w which act on fifth and eighth state
B_sig = [zeros(4,4);diag([1/Jm, 0, 0, 1/Ja3])];

C_sig = C;

% create systems for singular values
sys_ss_sig = ss(A_sig,B_sig,C_sig,[]);
sys_tf_sig = tf(sys_ss_sig);

% plot of the singular values, u to y
figure(2)
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
figure(3)
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

%% recreating the singular value plots from above manually with svd

% recreating plot of the singular values, u to y

sig_A1 = (sys_tf_sig(1,1));
[a_sysA1, b_sysA1, c_sysA1, d_sysA1] = ssdata(ss(sig_A1(1))); 
sig_A2 = ([(sys_tf_sig(1,1)) ; (sys_tf_sig(2,1)*s^2)]);
[a_sysA2, b_sysA2, c_sysA2, d_sysA2] = ssdata(ss(sig_A2)); 

w = logspace(-1,3,1e4);
Sig_A1 = zeros(length(w),1);
Sig_A2 = zeros(length(w),1);
for i = 1:1:length(w)
    out1 = c_sysA1*inv(1j*w(i)*eye(8)-a_sysA1)*b_sysA1;
    [u1,S1,v1] = svd(out1);
    out2 = c_sysA2*inv(1j*w(i)*eye(8)-a_sysA2)*b_sysA2;
    [u2,S2,v2] = svd(out2);
    Sig_A1(i) = S1(1,1);
    Sig_A2(i) = S2(1,1);
end

% frequency Response of Singular Values

figure(4)
loglog(w,Sig_A1,'b',w,Sig_A2,'r')
title('Singular Values from u to y, manual SVD')
xlabel('Frequency [rads/s]')
ylabel('Magnitude')
ylim([0.0001, 1000])
legend('y = q_m','y = [q_m Pdd]''')
grid on

% recreating plot of the singular values, w to y

sig_B1 = ([sys_tf_sig(1,1), sys_tf_sig(1,4)]);
[a_sysB1, b_sysB1, c_sysB1, d_sysB1] = ssdata(ss(sig_B1)); 
% sig_B2 = ([sys_tf_sig(1,1), sys_tf_sig(1,4) ;...
%     sys_tf_sig(2,1)*s^2,sys_tf_sig(2,4)*s^2]);
% [a_sysB2i, b_sysB2i, c_sysB2i, d_sysB2i] = ssdata(ss(sig_B2(2,2))); 

w = logspace(-1,3,1e4);
Sig_B1 = zeros(length(w),1);
Sig_B2 = zeros(length(w),1);
for i = 1:1:length(w)
    out1 = c_sysB1*inv(1j*w(i)*eye(8)-a_sysB1)*b_sysB1;
    [u1,S1,v1] = svd(out1);
%     out2 = c_sysB2i*inv(1j*w(i)*eye(8)-a_sysB2i)*b_sysB2i;
%     [u2,S2,v2] = svd(out2);
    Sig_B1(i) = S1(1,1);
%     Sig_B2(i) = S2(1,1);
end

% frequency Response of Singular Values

figure(5)
loglog(w,Sig_B1,'b'); % ,w,Sig_B2,'r')
title('Singular Values from w to y, manual SVD')
xlabel('Frequency [rads/s]')
ylabel('Magnitude')
ylim([0.0001, 1000])
legend('y = q_m'); % ,'y = [q_m Pdd]''')
grid on
hold off

%% References
%
% https://www.researchgate.net/publication/228715342_Robust_control_of_a_flexible_manipulator_arm_A_benchmark_problem
%