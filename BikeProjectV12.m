% Naji Husseini
% 9/26/2021
% BikeProjectV11.m
% Alec Johnson  Jose Rios  Daniella Sbordone  Kevin Zhang
% Code used to solve all sections of BMME 345 Project 1


clc
clear
close all

%% Declarations

% kinematic properties
th1 = -2*pi/3;
om2 = 1; %rad/s
al2 = 0; %rad/s^2
th2 = pi:pi/180:3*pi; %rad changing variable

% lengths
r1 = 30; % distance from hip to crank
r1 = r1*0.0254;  % m
r2 = 7; % pedal crank
r2 = r2*0.0254;  % m
r3 = 19.5; % lower leg
r3 = r3*0.0254;  % m
r4 = 19.75; % thigh
r4 = r4*0.0254;  % m

% mass calcs
tot_m = 156; %lbs
tot_m = tot_m/0.453592; % convert to kg
m_petal = 1; %kg, average weight of a pedal

m_thigh = tot_m*0.1;
m_LL = tot_m*0.061;


% Moment of Inertia
I_pc = (1/12)*m_petal*(r2)^2; % assume it is a rectangular bar
I_LL = m_LL*(r3*0.416)^2;
I_thigh = m_thigh*(r4*0.323)^2;

% link COMs
comL2 = 0.5;         % factor of COM of link 2 - petal crank

comL3_d = 0.394;         % distal factor of COM of link 3 - lower leg
comL3_p = 1 - comL3_d; % proximal factor of COM of link 3 

comL4_d = 0.567;         % distal factor of COM of link 4 - thigh
comL4_p = 1 - comL4_d; % proximal factor of COM of link 4 


T_hip_t = 125; %N*m  when crank is pi/2 to 3pi/2
T_hip_b = 0; %N*m when crank is at 3pi/2 to pi/2
g = -9.8; % Gravity, pretty important constant in kinetics section
% From this point onward, is my insertion


%% For Loop Calculations
% Insert this into main code file anywhere
%% Setup
kinval = zeros(6,361); % Initialize this to preallocate memory, to run faster
options = optimoptions('fsolve','Display','off');
guess = [2*pi/3 pi/4 0 0 0 0];
%% Solving
for k = 1:361
    sol = fsolve(@fourbar,guess,options,r1,r2,r3,r4,th1,th2(k),om2,al2);
    kinval(:,k) = sol;
    guess = sol;
end
%% Extracting values 
th3 = kinval(1,:); % Basically a matrix of stacked vectors, take each row
th4 = kinval(2,:); % to be a specific variable range. Minimizes for loop usage
om3 = kinval(3,:); % radians
om4 = kinval(4,:);
al3 = kinval(5,:);
al4 = kinval(6,:);

%% Deliverable #2: 
% Graphs theta 3 & 4 vs. theta 2, omega 3 & 4 vs. theta 2,
% and alpha 3 & 4 vs. theta 2

figure(1)
subplot(3, 1, 1)
title('\theta Values vs \theta_2')
hold on
plot(th2, th3)
plot(th2, th4)
xlabel('\theta_2 (radians)')
ylabel('\theta (rad)')
legend('\theta_3','\theta_4', 'Location', ' eastoutside')
hold off

subplot(3, 1, 2)
title('\omega Values vs \theta_2')
hold on
plot(th2,om3)
plot(th2,om4)
xlabel('\theta_2 (radians)')
ylabel('\omega (rad/sec)')
legend('\omega_3', '\omega_4', 'Location', ' eastoutside')
hold off

subplot(3, 1, 3)
title('\alpha Values vs \theta_2')
hold on
plot(th2, al3)
plot(th2, al4)
xlabel('\theta_3 (radians)')
ylabel('\alpha (rad/sec^2)')
legend('\alpha_3', '\alpha_4', 'Location', ' eastoutside')

xticks(pi:pi/4:3*pi);
xticklabels(["\pi" "5\pi/4" "3\pi/2" "7\pi/4" "2\pi" "9\pi/4" "5\pi/2" "11\pi/2" "3\pi"]);

hold off

% save to jpg
saveas(figure(1), 'kinematics.jpg')

%% Kinetics Setup
% Primary goal is to find the torque at the crank given the torque at the
% hip. Have to split it since the torque depends on which theta2 it is on.
% First, define position vectors from CoM to joints:

r12x = -comL2*r2*cos(th2); % A is defined as the joint at which the crank joins the bike frame
r12y = -comL2*r2*sin(th2); % Thus, assume it is a uniform bar for now.
r32x = -r12x; % This is simply the distance vector mirrored from the one
r32y = -r12y; % that goes from CoM to A.

r23x = -comL3_d*r3*cos(th3); % This is the foot + lower leg, and this vector
r23y = -comL3_d*r3*sin(th3); % points to the crank
r43x = comL3_p*r3*cos(th3); % Foot + lower leg, points towards knee
r43y = comL3_p*r3*sin(th3);

r34x = -comL4_d*r4*cos(th4); % Thigh, points towards knee
r34y = -comL4_d*r4*sin(th4);
r14x = comL4_p*r4*cos(th4); % Thigh, points to hip
r14y = comL4_p*r4*sin(th4);
% Also we need to determine linear accelerations so that they are defined
% at the heads of the vectors. We can scale them down later to find CoM
% linear accelerations. We use superposition to find the actual linear
% accelerations with respect to the frame.
a2x_alone = r2*(-al2.*sin(th2) - (om2).^2.*cos(th2));
a2y_alone = r2*(al2.*cos(th2) - (om2).^2 .*sin(th2));
a3x_alone = r3*(-al3.*sin(th3) - (om3).^2.*cos(th3));
a3y_alone = r3*(al3.*cos(th3) - (om3).^2.*sin(th3));
a4x_alone = r4*(-al4.*sin(th4) - (om4).^2.*cos(th4));
a4y_alone = r4*(al4.*cos(th4) - (om4).^2.*sin(th4));
% Below is the CoM calculations. Every variable here is already with
% respect to the frame and can be used in the matrix solve step.
a2x = a2x_alone*0.5;
a2y = a2y_alone*0.5;
a3x = a2x_alone + a3x_alone*comL3_d;
a3y = a2y_alone + a3y_alone*comL3_d;
a4x = a2x_alone + a3x_alone + a4x_alone*comL4_d;
a4y = a2y_alone + a3y_alone + a4y_alone*comL4_d;

%% Solving for Kinetics
% Remember that there are essentially three distinct regions since the thigh
% moment changes when the pedal is at the bottom or the top.
kin_ans = zeros(13,361); % Initialize matrix where answers will be to speed up
for k = 1:361 % Tp is the torque at the pedal
    %           F12x     F12y     F32x     F32y     Tp       F23x     F23y     F43x     F43y     F34x     F34y     F14x     F14y
    kinetics = [1        0        1        0        0        0        0        0        0        0        0        0        0   ;
                0        1        0        1        0        0        0        0        0        0        0        0        0   ;
                -r12y(k) r12x(k)  -r32y(k) r32x(k)  1        0        0        0        0        0        0        0        0   ;
                0        0        0        0        0        1        0        1        0        0        0        0        0   ;
                0        0        0        0        0        0        1        0        1        0        0        0        0   ;
                0        0        0        0        0        -r23y(k) r23x(k)  -r43y(k) r43x(k)  0        0        0        0   ;
                0        0        0        0        0        0        0        0        0        1        0        1        0   ;
                0        0        0        0        0        0        0        0        0        0        1        0        1   ;
                0        0        0        0        0        0        0        0        0        -r34y(k) r34x(k)  -r14y(k) r14x(k);
                0        0        1        0        0        1        0        0        0        0        0        0        0   ;
                0        0        0        1        0        0        1        0        0        0        0        0        0   ;
                0        0        0        0        0        0        0        1        0        1        0        0        0   ;
                0        0        0        0        0        0        0        0        1        0        1        0        0   ];
    values = [m_petal*a2x(k) m_petal*(a2y(k)-g) I_pc*al2 ...
              m_LL*a3x(k) m_LL*(a3y(k)-g) I_LL*al3(k) ...
              m_thigh*a4x(k) m_thigh*(a4y(k)-g) I_thigh*al4(k)-T_hip_t ...
              0 0 0 0]'; % Last four 0s represent Newton's laws
    kin_ans(:,k) = kinetics\values;
end
for k = 92:271 % Overwrites middle portion of kin_ans such that it accounts for
               % the hip torque being zero from 270 degrees to 430 degrees
    kinetics = [1        0        1        0        0        0        0        0        0        0        0        0        0   ;
                0        1        0        1        0        0        0        0        0        0        0        0        0   ;
                -r12y(k) r12x(k)  -r32y(k) r32x(k)  1        0        0        0        0        0        0        0        0   ;
                0        0        0        0        0        1        0        1        0        0        0        0        0   ;
                0        0        0        0        0        0        1        0        1        0        0        0        0   ;
                0        0        0        0        0        -r23y(k) r23x(k)  -r43y(k) r43x(k)  0        0        0        0   ;
                0        0        0        0        0        0        0        0        0        1        0        1        0   ;
                0        0        0        0        0        0        0        0        0        0        1        0        1   ;
                0        0        0        0        0        0        0        0        0        -r34y(k) r34x(k)  -r14y(k) r14x(k);
                0        0        1        0        0        1        0        0        0        0        0        0        0   ;
                0        0        0        1        0        0        1        0        0        0        0        0        0   ;
                0        0        0        0        0        0        0        1        0        1        0        0        0   ;
                0        0        0        0        0        0        0        0        1        0        1        0        0   ];
    values = [m_petal*a2x(k) m_petal*(a2y(k)-g) I_pc*al2 ...
              m_LL*a3x(k) m_LL*(a3y(k)-g) I_LL*al3(k) ...
              m_thigh*a4x(k) m_thigh*(a4y(k)-g) I_thigh*al4(k) ...
              0 0 0 0]';
    kin_ans(:,k) = kinetics\values;
end

%% Deliverable #3:
% Graphs torque at the hip and torque at the pedal vs. theta 2
% throughout one rotation of the pedal about the frame from pi
% to 3*pi

figure(2)
hold on;
th2_graph125 = th2(92:271);
th2_graph01 = th2(1:91);
th2_graph02 = th2(272:end);
T_hip_0 = zeros(size(th2_graph125));
T_hip_01 = ones(size(th2_graph01))*T_hip_t;
T_hip_02 = ones(size(th2_graph02))*T_hip_t;

T_hip_final = reshape([T_hip_01 T_hip_0 T_hip_02].',[],1);

plot(th2,kin_ans(5,:), th2, T_hip_final ); 
xlabel('\theta_2 (rad)')
ylabel('Torque at the pedal (N*m)')
title('Torque at the Hip & Pedal vs. Theta 2')
legend('Pedal Torque','Hip Torque')
xticks(pi:pi/4:3*pi);
xticklabels(["\pi" "5\pi/4" "3\pi/2" "7\pi/4" "2\pi" "9\pi/4" "5\pi/2" "11\pi/2" "3\pi"]);
hold off;

% save to jpg
saveas(figure(2), 'torque.jpg')


%% Components of Force
Ftx = kin_ans(10,:); % x-components of F34x
Fty = kin_ans(11,:); % y-components of F34y
transform = zeros(2,361); % initialize matrix
th_4 = kinval(2,:); % grab the theta 4 values
for k = 1:361
    T = [cos(th_4(k)) sin(th_4(k)); -sin(th_4(k)) cos(th_4(k))]; % Transformation Matrix
    calc = T*[Ftx(k);Fty(k)]; % Transform vectors into rotating coordinate system
    transform(:,k) = calc; % Drop values in matrix
    % Row one is x components
    % Row two is y components
end

%% Deliverable #4:
% Graph of force components at knee in the direction of the 
% thigh and perpendicular to the thigh vs. theta 2 

figure(3)
plot(th2,transform(1,:),'r', th2, transform(2,:),'b'); 
xlabel('\theta_2 (rad)')
ylabel('Force at the Knee (N)')
title('Force at Knee vs. \theta_2')
legend('Force in Thigh Direction','Force Perpendicular to Thigh')
xticks(pi:pi/4:3*pi);
xticklabels(["\pi" "5\pi/4" "3\pi/2" "7\pi/4" "2\pi" "9\pi/4" "5\pi/2" "11\pi/2" "3\pi"]);

% save to jpg
saveas(figure(3), 'kneeForces.jpg')


% %% Deliverable #5:
% % Stick‐figure video showing the links and a trace 
% % of the knee’s position
% % Preallocation:
% mem = zeros(1,361);
% pxA = mem;
% pyA = mem;
% pxB = mem;
% pyB = mem;
% pxC = mem;
% pyC = mem;
% pxD = mem;
% pyD = mem;
% % Position equations
% for k = 1:361
%     pxA(k) = r1*cos(th1); % position of link b/w frame and pedal crank
%     pyA(k) = r1*sin(th1);
%     
%     pxB(k) = r1*cos(th1) + r2*cos(th2(k)); % position of pedal
%     pyB(k) = r1*sin(th1) + r2*sin(th2(k));
%     
%     pxC(k) = r1*cos(th1) + r2*cos(th2(k)) + r3*cos(th3(k)); % position of knee
%     pyC(k) = r1*sin(th1) + r2*sin(th2(k)) + r3*sin(th3(k));
%     
%     pxD(k) = r1*cos(th1) + r2*cos(th2(k)) + r3*cos(th3(k)) + r4*cos(th4(k));  % position of hip
%     pyD(k) = r1*sin(th1) + r2*sin(th2(k)) + r3*sin(th3(k)) + r4*sin(th4(k));
%     
% end
% 
% % initialize video
% myVideo = VideoWriter('cycling'); %open video file
% myVideo.FrameRate = 30;
% open(myVideo)
% 
% % Creates video figure
% for k = 1:361
%     figure(4)
%     plot(pxC, pyC, 'g','MarkerSize', 30) % plots knee range
%     hold on
%     plot([pxA(k) pxB(k)], [pyA(k) pyB(k)], 'b','MarkerSize', 30) % plots links
%     plot([pxB(k) pxC(k)], [pyB(k) pyC(k)], 'r','MarkerSize', 30)
%     plot([pxC(k) pxD(k)], [pyC(k) pyD(k)], 'MarkerSize', 30)
%     plot([pxD(k) pxA(k)], [pyD(k) pyA(k)], 'MarkerSize', 30)
%     title('Pedaling System')
%     ylabel('Vertical Position (m)')
%     xlabel('Horizontal Position (m)')
%     hold off
%     axis([-1 0.5 -1 0.5])
%     legend('knee', 'pedal', 'lower leg', 'thigh', 'frame')
%     frame = getframe(gcf); %get frame
%     writeVideo(myVideo, frame);
% end
% close(myVideo)