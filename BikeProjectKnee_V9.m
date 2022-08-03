% BikeProjectKnee_V8.m
% 11/17/2021
% Alec Johnson  Jose Rios  Daniella Sbordone  Kevin Zhang
% Code used to solve all sections of BMME 345 Project 2


clc
clear
close all

%% Declarations

% Global variables
n = 3; % Number of cycles
f32y = csvread('f32y.csv'); % read in f32y
t = f32y(:,1)';
F32y = f32y(:,2)';
% Bike dimensions
% lengths
in2m = 0.0254;
seatHeight = 15*in2m; % height from crank to seat (m)
seatLength = 17*in2m; % horizontal distance from crank to seat (m)
m2 = 0.5; % pedal mass (kg)
h2 = 1*in2m; % pedal height (m)

% Dimensions
r1 = sqrt(seatHeight^2+seatLength^2);
r2 = 4*in2m;
r3 = 21*in2m;
r4 = 19*in2m;

% Link 1 (frame) kinematics
th1 = -pi/2 - atan(seatLength/seatHeight);
om1 = 0;
al1 = 0;
% om2 defined later
al2 = 0; %rad/s^2
th2 = linspace(pi,(2*n+1)*pi,length(F32y)); %rad changing variable

% mass calcs
tot_m = 75; %kg
m_petal = 1; %kg, average weight of a pedal

m_thigh = tot_m*0.1;
m_LL = tot_m*0.061; % Includes foot and lower leg


% Moment of Inertia (kg*m^2)
I_pc = (1/12)*m2*(h2^2 + r2^2); % assume it is a rectangular bar
I_LL = m_LL*(r3*0.416)^2;
I_thigh = m_thigh*(r4*0.323)^2;

% link COMs (m)
comL2 = 0.5;         % factor of COM of link 2 - petal crank

comL3_d = 0.394;         % distal factor of COM of link 3 - lower leg
comL3_p = 1 - comL3_d; % proximal factor of COM of link 3 

comL4_d = 0.567;         % distal factor of COM of link 4 - thigh
comL4_p = 1 - comL4_d; % proximal factor of COM of link 4 

g = -9.81; % Gravity, pretty important constant in kinetics section

% Implant values
implant_len = 75*1e-3;

%% Finding om2
F_avg = mean(F32y);
F_bool_low = F32y <= F_avg + 2; % Thresholding method
F_bool_high = F32y >= F_avg - 2; 
F_bool = and(F_bool_low,F_bool_high); % Find time values where force close to average
t_points = t(F_bool); % indexing
period = 2*mean(diff(t_points)); % Find average half periods, x2 to get period
om2 = 2*pi/period; % angular frequency found
%% Setup
elements = 1:length(F32y);
kinval = zeros(6,length(F32y)); % Initialize this to preallocate memory, to run faster
options = optimoptions('fsolve','Display','off');
guess = [2*pi/3 pi/4 0 0 0 0];
%% Solving
for k = elements
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


%% Kinetics Setup
%Using F32, solve for both the torque at the pedal and hip

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
kin_ans = zeros(13,length(elements)); % Initialize matrix where answers will be to speed up
for k = elements
    %           F12x     F12y     F32x     Thip     Tp       F23x     F23y     F43x     F43y     F34x     F34y     F14x     F14y
    kinetics = [1        0        1        0        0        0        0        0        0        0        0        0        0   ;
                0        1        0        0        0        0        0        0        0        0        0        0        0   ;
                -r12y(k) r12x(k)  -r32y(k) 0        1        0        0        0        0        0        0        0        0   ;
                0        0        0        0        0        1        0        1        0        0        0        0        0   ;
                0        0        0        0        0        0        1        0        1        0        0        0        0   ;
                0        0        0        0        0        -r23y(k) r23x(k)  -r43y(k) r43x(k)  0        0        0        0   ;
                0        0        0        0        0        0        0        0        0        1        0        1        0   ;
                0        0        0        0        0        0        0        0        0        0        1        0        1   ;
                0        0        0        1        0        0        0        0        0        -r34y(k) r34x(k)  -r14y(k) r14x(k);
                0        0        1        0        0        1        0        0        0        0        0        0        0   ;
                0        0        0        0        0        0        1        0        0        0        0        0        0   ;
                0        0        0        0        0        0        0        1        0        1        0        0        0   ;
                0        0        0        0        0        0        0        0        1        0        1        0        0   ];
    values = [m_petal*a2x(k) m_petal*(a2y(k)-g)-F32y(k) I_pc*al2-r32x(k)*F32y(k) ...
              m_LL*a3x(k) m_LL*(a3y(k)-g) I_LL*al3(k) ...
              m_thigh*a4x(k) m_thigh*(a4y(k)-g) I_thigh*al4(k) ...
              0 -F32y(k) 0 0]';
    kin_ans(:,k) = kinetics\values;
end

%% Deliverable #2

figure(2)
plot(th2,kin_ans(4,:));
hold on;
plot(th2,kin_ans(5,:));
legend('Thip','Tpedal');
xlabel('\theta_2 (radians)');
ylabel('Torque (N*m)');
title('Torque vs. \theta_2');
xticks(pi:pi/2:7*pi);
xticklabels(["\pi" "3\pi/2" "2\pi" "5\pi/2" "3\pi" "7\pi/2" ...
    "4\pi" "9\pi/2" "5\pi" "11\pi/2" "6\pi" "13\pi/2" "7\pi"]);
hold off;

% save to jpg
saveas(figure(2), 'torque.jpg')

%% Components of Force


Ftx = kin_ans(8,:); % x-components of F34x
Fty = kin_ans(9,:); % y-components of F34y
transform = zeros(2,length(elements)); % initialize matrix
th_3 = kinval(1,:); % grab the theta 4 values
for k = elements
    T = [cos(th_3(k)) sin(th_3(k)); -sin(th_3(k)) cos(th_3(k))]; % Transformation Matrix
    calc = T*[Ftx(k);Fty(k)]; % Transform vectors into rotating coordinate system in frame with all axes
    transform(:,k) = calc; % Drop values in matrix
    % Row one is x components
    % Row two is y components
end

%% Deliverable #3:
% Graph of force components at knee in the direction of the 
% lower leg and perpendicular to the lower leg vs. theta 2 

figure(3)
plot(th2,transform(1,:),'r', th2, transform(2,:),'b'); 
xlabel('\theta_2 (rad)')
ylabel('Force at the Knee (N)')
title('Force at Knee vs. \theta_2')
legend('Force in Lower Leg Direction','Force Perpendicular to Lower Leg','Location', 'southoutside')
xticks(pi:pi/2:7*pi);
xticklabels(["\pi" "3\pi/2" "2\pi" "5\pi/2" "3\pi" "7\pi/2" ...
    "4\pi" "9\pi/2" "5\pi" "11\pi/2" "6\pi" "13\pi/2" "7\pi"]);

% save to jpg
saveas(figure(3), 'kneeForces.jpg')

%% Declarations for Stress
% Define important constants like E, Poisson's ratio, etc. The force
% parallel the knee was rotated until it was parallel with the x-axis.

E = 113.8e9; % in Pa, work with all SI units
v = 0.342;
R_stem = 6.5e-3; % 6.5 mm in meters, corresponds to 13 mm diameter Depuy PFC Sigma
xarea_stem = pi*R_stem^2; % Assume this is a cylinder
stress_x = transform(1,:)/xarea_stem; % Normal stress
stress_x_MPa = stress_x / 1e6;
shear_xy = transform(2,:)/xarea_stem; % Shear stress
shear_xy_MPa = shear_xy / 1e6;
sigma_3_values = zeros(1,length(elements)); % Maximum compressive stress
for k = elements
    sigma = [stress_x(k) shear_xy(k);
             shear_xy(k)     0      ];
    sigma_3_values(k) = min(eig(sigma),[],'all');
end
sigma_3_MPa = sigma_3_values / 1e6;
%% Deliverable #4
%Plots for normal, shear, and max compressive principle stresses in the 
%stem of the knee implant

figure(4)
subplot(3, 1, 1)
title('Normal \sigma vs \theta_2 in the Stem')
hold on
plot(th2, stress_x_MPa,'r')
xlabel('\theta_2 (radians)')
ylabel('\sigma_N (MPa)')
xticks(pi:pi/2:7*pi);
xticklabels(["\pi" "3\pi/2" "2\pi" "5\pi/2" "3\pi" "7\pi/2" ...
    "4\pi" "9\pi/2" "5\pi" "11\pi/2" "6\pi" "13\pi/2" "7\pi"]);
hold off
subplot(3, 1, 2)
title('Shear Stress vs \theta_2 in the Stem')
hold on
plot(th2,shear_xy_MPa,'m')
xlabel('\theta_2 (radians)')
ylabel('\tau_x_y (MPa)')
xticks(pi:pi/2:7*pi);
xticklabels(["\pi" "3\pi/2" "2\pi" "5\pi/2" "3\pi" "7\pi/2" ...
    "4\pi" "9\pi/2" "5\pi" "11\pi/2" "6\pi" "13\pi/2" "7\pi"]);
hold off

subplot(3, 1, 3)
title('Compressive \sigma_p vs \theta_2 in the Stem')
hold on
plot(th2,sigma_3_MPa,'b')
xlabel('\theta_3 (radians)')
ylabel('\sigma_P (MPa)')
xticks(pi:pi/2:7*pi);
xticklabels(["\pi" "3\pi/2" "2\pi" "5\pi/2" "3\pi" "7\pi/2" ...
    "4\pi" "9\pi/2" "5\pi" "11\pi/2" "6\pi" "13\pi/2" "7\pi"]);

% save to jpg
saveas(figure(4), 'stress.jpg')

%% Declarations for Strain

G = E/(2*(1+v)); % Shear modulus defined
strain_x = stress_x / E; % Convert stress to strain
strain_y = -v*strain_x;
shearstrain_xy = shear_xy / G;
shearstrain_elastic = shearstrain_xy/2; % Converting to elastic strain for calc
tau_3 = zeros(1,length(elements));
for k = elements
    ep = [strain_x(k) shearstrain_elastic(k);
          shearstrain_elastic(k) strain_y(k)];
    tau_3(k) = min(eig(ep),[],'all');
end


%% Deliverable #5
%Plots for normal, shear, and max compressive principle strains in the 
%stem of the knee implant

figure(5)
subplot(3, 1, 1)
title('Normal \epsilon vs \theta_2 in the Stem')
hold on
plot(th2, 1e6*strain_x,'r') % Plot in microstrain
plot(th2, 1e6*strain_y);
xlabel('\theta_2 (radians)')
ylabel('\epsilon_N (\mu\epsilon)') 
xticks(pi:pi/2:7*pi);
xticklabels(["\pi" "3\pi/2" "2\pi" "5\pi/2" "3\pi" "7\pi/2" ...
    "4\pi" "9\pi/2" "5\pi" "11\pi/2" "6\pi" "13\pi/2" "7\pi"]);hold off

subplot(3, 1, 2)
title('Shear Strain vs \theta_2 in the Stem')
hold on
plot(th2,1e6*shearstrain_xy,'m') % Plot in microstrain
xlabel('\theta_2 (radians)')
ylabel('\gamma_x_y (\mu\epsilon)') 
xticks(pi:pi/2:7*pi);
xticklabels(["\pi" "3\pi/2" "2\pi" "5\pi/2" "3\pi" "7\pi/2" ...
    "4\pi" "9\pi/2" "5\pi" "11\pi/2" "6\pi" "13\pi/2" "7\pi"]);
hold off

subplot(3, 1, 3)
title('Compressive \epsilon_p vs \theta_2 in the Stem')
hold on
plot(th2,1e6*tau_3,'b') % Plot in microstrain
xlabel('\theta_3 (radians)')
ylabel('\epsilon_P (\mu\epsilon)')
xticks(pi:pi/2:7*pi);
xticklabels(["\pi" "3\pi/2" "2\pi" "5\pi/2" "3\pi" "7\pi/2" ...
    "4\pi" "9\pi/2" "5\pi" "11\pi/2" "6\pi" "13\pi/2" "7\pi"]);

% save to jpg
saveas(figure(5), 'strain.jpg')
