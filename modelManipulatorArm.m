%% ECSE 6460 Multivariable Control
% Final Project - Robust Control of a Flexible Manipulator
%
% Model
%
% Last Edited: April 10, 2017

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

% Note for linearized model, k1 = constant and Td = 0
% k1 is the nonlinear stiffness for the gear box

Jm = 5e-3; Ja1 = 2e-3; Ja2 = 0.02; Ja3 = 0.02; % moment of intertia

k1_low = 16.7;
k1_high = 100;
k1 = (100+16.7)/2;
k1_nom = linspace(k1_low, k1_high, 5);

k2 = 110; k3 = 80; % stiffness
d1 = 0.08; d2 = 0.06; d3 = 0.08; % damping
fm = 6e-3; fa1 = 1e-3; fa2 = 1e-3;fa3 = 1e-3; % viscous friction
n = 220; % gear ratio
l1 = 20e-3; l2 = 600e-3; l3 = 1530e-3; % link lengths
Td = 0; %0.5e-3; % time delay


%% Linearized Model

J = diag([Jm, Ja1, Ja2, Ja3]);

D = [d1, -d1, 0, 0;
    -d1, (d1+d2), -d2, 0;
    0, -d2, (d2+d3), -d2;
    0, 0, -d3, d3];

F = diag([fm, fa1, fa2, fa3]);

K = [k1, -k1, 0, 0;
    -k1, (k1+k2), -k2, 0;
    0, -k2, (k2+k3), -k3;
    0, 0, -k3, k3];

A = [zeros(4,4), eye(4);
    -inv(J)*K, -inv(J)*(D+F)];

B = [zeros(4,4);
    inv(J)];

E = (1/n) * [0, l1, l2, l3, zeros(1,4)];
C = [1, zeros(1,7);
    E];

sys_ss = ss(A,B,C,[]);

%% Transfer function, motor torque to motor acceleration, nominal moderl
%
% 

figure(1);
grid on
hold on
for i = 1:length(k1_nom)
    
    K_nom = [k1_nom(i), -k1_nom(i), 0, 0;
        -k1_nom(i), (k1_nom(i)+k2), -k2, 0;
        0, -k2, (k2+k3), -k3;
        0, 0, -k3, k3];
    A_nom = [zeros(4,4), eye(4);
        -inv(J)*K_nom, -inv(J)*(D+F)];
    
    C_nom = [0,0,0,0,1,0,0,0];
    
    
    sys_ss_nom = ss(A_nom,B,C_nom,[]);
    
    sys_tf_nom = tf(sys_ss_nom(1,1))*s;
    
    h = bodeplot(sys_tf_nom);
    setoptions(h,'FreqUnits','Hz','PhaseVisible','off');
    xlim([1, 100])
    
end
title('Transfer Function, Motor Torque to Motor Acceleration')
legend(['k_{1} = ',num2str(k1_nom(1))],['k_{1} = ',num2str(k1_nom(2))],['k_{1} = ',num2str(k1_nom(3))],['k_{1} = ',num2str(k1_nom(4))],['k_{1} = ',num2str(k1_nom(5))])
hold off


%% Investigate Singular Values - Using K1 = max = 100

% Hinf Norm from the Singular Value Decomposition

% using Skogestad Eq. 3.49, pp. 81
% Hinf norm is defined as the peak of the maximum singular value of the
% frequency response.

K_sin = [k1_nom(5), -k1_nom(5), 0, 0;
    -k1_nom(5), (k1_nom(5)+k2), -k2, 0;
    0, -k2, (k2+k3), -k3;
    0, 0, -k3, k3];

A_sin = [zeros(4,4), eye(4);
    -inv(J)*K_sin, -inv(J)*(D+F)]

B_sin_u2Yq = [zeros(4,4);diag([1/Jm, 0, 0, 0])]; % u to y
B_sin_w2Yq = [zeros(4,4);diag([1/Jm, 0, 0, 1/Ja3])]; % w to y


C_sin_q = [1,0,0,0,0,0,0,0]; % y = qm (blue)
C_sin_P = [1, zeros(1,7); E]; % y = [qm Pdd]' (red)

% singular values from u to y
sys_ss_sin_u2Yqq = ss(A_sin,B_sin_u2Yq,C_sin_q,[]); % y = qm (blue) 
sys_ss_sin_u2YqP = ss(A_sin,B_sin_u2Yq,C_sin_P,[]);
sys_tf_sin_u2YqP = tf(sys_ss_sin_u2YqP)

figure(2);
sigmaplot(sys_ss_sin_u2Yqq,{1e-1,1e4},'b')
title('Singular Values from u to y')
xlabel('Frequency [rads/s]')
ylabel('Magnitude [dB]')
grid on

sys_ss_sin_w2Yqq = ss(A_sin,B_sin_w2Yq,C_sin_q,[]);
figure(3);
sigmaplot(sys_ss_sin_w2Yqq,{1e-1,1e4},'b')
title('Singular Values from w to y')
xlabel('Frequency [rads/s]')
ylabel('Magnitude [dB]')
grid on

% w = logspace(-1,3,1e4);
% s1 = zeros(length(w),1);
% s2 = zeros(length(w),1);
% for i = 1:1:length(w)
%     out = C_sin_U2Y*inv(1j*w(i)*eye(8)-A_sin)*B_sin_U2Y;
%     [u,s,v] = svd(out)
%     s1(i) = s(1,1);
%     %s2(i) = s(2,2);
% end
% % 
% Hinf_norm = max(s1(:))
% 
% % frequency Response of Singular Values
% 
% % % % figure(2)
% % % % loglog(w,s1,'b');%,w,s2,'r')
% % % title('Frequency Response of Singular Values')
% % % xlabel('Frequency [rads/s]')
% % % ylabel('Magnitude')
% % % grid on
% % % %legend('\sigma_{1}', '\sigma_{2}')
% % % 
% % % figure(3);
% % % sigmaplot(sys_ss_sin_U2Y,{1e-1,1e4})











% % %% Transfer function Matrix
% % 
% % sys_tf = tf(sys_ss)
% % 
% % %% Eigenvalues and Eigenvectors
% % 
% % [V,D,W] = eig(A);
% % 
% % fprintf('EigenValues: \n\n');
% % disp(D*ones(8,1));
% % fprintf('Right EigenVectors: \n\n');
% % disp(V);
% % fprintf('Left EigenVectors: \n\n');
% % disp(W);
% % 
% % %% Input / Output Pole Vectors
% % 
% % % all poles
% % all_poles = pole(sys_ss)
% % 
% % % using Skogestad Table 4.2: Matlab Commands to Find Pole Vectors
% % [Q,Pi] = eig(A');
% % UP = B' * Q; % input pole vectors
% % [T,Po] = eig(A);
% % YP = C * T; % output pole vectors
% % Shouldbezero = Po-Pi % if not, pole vectors refer to different poles
% % 
% % fprintf('input pole vectors: UP_transpose \n\n')
% % disp(UP')
% % fprintf('output pole vectors: YP_transpose \n\n')
% % disp(YP')
% % 
% % %% Zeros and Zero Directions
% % 
% % % invariant zeros, MATLAB Docs recommend using ss model
% % fprintf('All Zeros \n')
% % z0 = tzero(sys_ss)
% % disp('No system zeros. Therefore no zero directions.')
% % 
% % %% Controllability  and Observability
% % 
% % % state controllability
% % 
% % % using row rank of controllability matrix
% % Co = ctrb(A,B);
% % fprintf('\nRank of Controllability Matrix')
% % rCo = rank(Co) % row rank
% % unCo = length(A) - rCo % number of uncontrollable states
% % 
% % % using Grammian
% % fprintf('\nRank of Controllability Grammian')
% % Wc = gram(sys_ss,'c') 
% % rWc = rank(Wc)
% % unWc = length(A) - rWc
% % 
% % % Popov-Belevitch-Hautus (PBH) Test
% % lambda_x = eig(A);
% % unPBH_c = ones(length(lambda_x),1);
% % for i = 1:1:length(lambda_x)
% %     PBH_c = [A-lambda_x(i)*eye(8), B];
% %     rPBH_c = rank(PBH_c);
% %     unPBH_c(i) = length(A) - rPBH_c;
% % end
% % clear i
% % fprintf('\nPBH Test of Controllability performed in script \n')
% % 
% % if unCo ~= 0 | unWc ~= 0 | unPBH_c ~= 0
% %     fprintf('\n!!!!!Uncontrollable Flag!!!!! \n')
% % end
% % 
% % % state observability
% % 
% % % using rank of kernal of observability matrix.
% % Ob = obsv(A,C);
% % rOb = rank(Ob) % column rank
% % unOb = length(A) - rOb % number of unobservable states
% % 
% % % % using Grammian
% % % fprintf('\nRank of Observability Grammian')
% % % Wo = gram(sys_ss,'o') 
% % % rWo = rank(Wo)
% % % unWo = length(A) - rWo
% % 
% % % Popov-Belevitch-Hautus (PBH) Test
% % lambda_y = eig(A);
% % unPBH_o = ones(length(lambda_y),1);
% % for i = 1:1:length(lambda_y)
% %     PBH_o = [A-lambda_y(i)*eye(8), B];
% %     rPBH_o = rank(PBH_o);
% %     unPBH_o(i) = length(A) - rPBH_o;
% % end
% % clear i
% % fprintf('\nPBH Test of Observability performed in script')
% % 
% % % if unOb ~= 0 | unWo ~= 0 | unPBH_o ~= 0
% % %     disp('UnObservable Flag')
% % % end
% % % fprintf('\nIf no flag raised, Observable\n')
% % 
% % 
% % %% H2 Norm from the State Space Representation
% % 
% % % using Zhou Eq. 4.1, pp. 53
% % H2_norm  = sqrt(trace(C*Wc*C')) 
% % 
% % %% Hinf Norm from the Singular Value Decomposition
% % 
% % % using Skogestad Eq. 3.49, pp. 81
% % % Hinf norm is defined as the peak of the maximum singular value of the
% % % frequency response.
% % 
% % w = logspace(-1,3,1e4);
% % s1 = zeros(length(w),1);
% % s2 = zeros(length(w),1);
% % for i = 1:1:length(w)
% %     out = C*inv(1j*w(i)*eye(8)-A)*B;
% %     [u,s,v] = svd(out);
% %     s1(i) = s(1,1);
% %     s2(i) = s(2,2);
% % end
% % 
% % Hinf_norm = max(s1(:))
% % 
% % %% frequency Response of Singular Values
% % 
% % figure(1)
% % loglog(w,s1,'b',w,s2,'r')
% % title('Problem 1-8: Frequency Response of Singular Values')
% % xlabel('Frequency [rads/s]')
% % ylabel('Magnitude')
% % grid on
% % legend('\sigma_{1}', '\sigma_{2}')

%% References
%
% https://www.researchgate.net/publication/228715342_Robust_control_of_a_flexible_manipulator_arm_A_benchmark_problem
%