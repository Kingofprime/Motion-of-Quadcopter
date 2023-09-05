clc; clear; clear all;
%% Inital Conditions

g= 9.81;
m=0.450;
l=0.225;
k=2.98e-7;
b = 1.14e-6;
ixx=4.85e-3;
iyy=ixx;
izz= 8.8e-3;
I=[ixx,0,0; 0, iyy,0; 0,0, izz];
r0vec=[0,0,1]';
rdot0=[0,0,0]';
psi0=deg2rad(10);
theta0=deg2rad(10);
phi0=deg2rad(10);
p0=0;
q0=0;
r0=0;
zr=10;
phi_r=0;
theta_r=0;
psi_r=0;
zdot_r=0;
phidot_r=0;
thetadot_r=0;
psidot_r=0;
x0=[r0vec; rdot0; psi0; theta0; phi0; p0; q0; r0;zr; phi_r; theta_r; psi_r;zdot_r; phidot_r;thetadot_r;psidot_r];
t=120;
tspan = linspace(0,t,10001);

%% set ode45

options = odeset('AbsTol',1e-12,'RelTol',1e-12);
[t2,p2]=ode45(@(t,x)ode(t,x,g,m,l,k,b,I,zr,phi_r,theta_r,psi_r, zdot_r,phidot_r,thetadot_r,psidot_r),tspan,x0,options);

%% Inertial Frame to BODY Frame

for i=1:length(t2)
    dcm1= [1,0,0 ; 0, cos(p2(i,7)), -sin(p2(i,7)); 0 , sin(p2(i,7)), cos(p2(i,7))];
    dcm2=[cos(p2(i,8)) , 0 , sin(p2(i,8)); 0,1,0; -sin(p2(i,8)), 0 , cos(p2(i,8))];
    dcm3 = [cos(p2(i,9)) , -sin(p2(i,9)), 0; sin(p2(i,9)), cos(p2(i,9)), 0; 0,0,1];
    dcm=dcm1*dcm2*dcm3;
    bf(i,:)=dcm*[p2(i,4);p2(i,5);p2(i,6)];
end

%% rotor angular velocities versus time

for i=1:length(t2)
    T = (g+(zdot_r-p2(i,6))+(zr-p2(i,3)))*m/(cos(p2(i,7))*cos(p2(i,8)));
    A= [1,0, -sin(p2(i,8)) ; 0, cos(p2(i,7)), cos(p2(i,8))*sin(p2(i,7)); 0, -sin(p2(i,7)), cos(p2(i,8))*cos(p2(i,7))];
    eulerrates=A\(p2(i,10:12))';
    L = I(1,1)*((phidot_r-eulerrates(1))+(phi_r-p2(i,7)));
    M = I(2,2)*((thetadot_r-eulerrates(2))+(theta_r-p2(i,8)));
    N = I(3,3)*((psidot_r-eulerrates(3))+(psi_r-p2(i,9)));

    om1=sqrt((T/(4*k))-(M/(2*k*l))+(N/(4*b)));
    om2=sqrt((T/(4*k))-(L/(2*k*l))-(N/(4*b)));
    om3=sqrt((T/(4*k))+(M/(2*k*l))+(N/(4*b)));
    om4=sqrt((T/(4*k))+(L/(2*k*l))-(N/(4*b)));

    om(i,:)=[om1;om2;om3;om4];
end
%% plots
% Position vs Time
figure
subplot(3,1,1)
plot(t2,p2(:,1));
ylabel('X position (m)');
xlabel('time (sec)');
title('X position in intertial frame vs time');
subplot(3,1,2)
plot(t2,p2(:,2));
ylabel('Y position (m)');
xlabel('time (sec)');
title('Y position in intertial frame vs time');
subplot(3,1,3)
plot(t2,p2(:,3));
ylabel('Z position (m)');
xlabel('time (sec)');
title('Z position in intertial frame vs time');

% Velocity vs Time
figure
subplot(3,1,1)
plot(t2,bf(:,1));
ylabel('X velocity (m\s)');
xlabel('time (sec)');
title('X velocity in body frame vs time');
subplot(3,1,2)
plot(t2,bf(:,2));
ylabel('Y velocity (m\s)');
xlabel('time (sec)');
title('Y velocity in body frame vs time');
subplot(3,1,3)
plot(t2,bf(:,3));
ylabel('Z velocity (m\s)');
xlabel('time (sec)');
title('Z velocity in body frame vs time');

% Euler Angle vs Time
figure
subplot(3,1,1)
plot(t2,p2(:,7));
ylabel('phi (rad)');
xlabel('time (sec)');
title('Euler angle (phi) vs time');
subplot(3,1,2)
plot(t2,p2(:,8));
ylabel('theta (rad)');
xlabel('time (sec)');
title('Euler angle (theta) vs time');
subplot(3,1,3)
plot(t2,p2(:,9));
ylabel('psi (rad)');
xlabel('time (sec)');
title('Euler angle (psi) vs time');

% Angular Velocity vs Time
figure
subplot(3,1,1)
plot(t2,p2(:,10));
ylabel('p (rad/sec)');
xlabel('time');
title('Angular Velocity (p) vs time (sec)');
subplot(3,1,2)
plot(t2,p2(:,11));
ylabel('q (rad/sec)');
xlabel('time (sec)');
title('Angular Velocity (q) vs time');
subplot(3,1,3)
plot(t2,p2(:,12));
ylabel('r (rad/sec)');
xlabel('time (sec)');
title('Angular Velocity (r) vs time');

% Autopilot omega vs time
figure
subplot(4,1,1)
plot(t2,om(:,1));
ylabel('om1 (rad/sec)');
xlabel('time (sec)');
title('Auto Pilot omega 1 vs time');
subplot(4,1,2)
plot(t2,om(:,2));
ylabel('om2 (rad/sec)');
xlabel('time (sec)');
title('Auto Pilot omega 2 vs time');
subplot(4,1,3)
plot(t2,om(:,3));
ylabel('om3 (rad/sec)');
xlabel('time (sec)');
title('Auto Pilot omega 3 vs time');
subplot(4,1,4)
plot(t2,om(:,4));
ylabel('om4 (rad/sec)');
xlabel('time (sec)');
title('Auto Pilot omega 4 vs time');


%% ODE function

function [xdot] = ode(t,x,g,m,l,k,b,I,zr,phi_r,theta_r,psi_r, zdot_r,phidot_r,thetadot_r,psidot_r)

T = (g+(zdot_r-x(6))+(zr-x(3)))*m/(cos(x(7))*cos(x(8)));
A= [1,0, -sin(x(8)) ; 0, cos(x(7)), cos(x(8))*sin(x(7)); 0, -sin(x(7)), cos(x(8))*cos(x(7))];
eulerrates=A\x(10:12);
L = I(1,1)*((phidot_r-eulerrates(1))+(phi_r-x(7)));
M = I(2,2)*((thetadot_r-eulerrates(2))+(theta_r-x(8)));
N = I(3,3)*((psidot_r-eulerrates(3))+(psi_r-x(9)));

omega = x(10:12);
tau=[L;M;N];
omegadot = I\(cross(-omega,I*omega)+tau);

% output
xdot=zeros(size(x));
xdot(1:3)=x(4:6);
xdot(4)= T/m*(cos(x(9))*sin(x(8))*cos(x(7))+sin(x(9))*sin(x(7)));
xdot(5)= T/m*(sin(x(9))*sin(x(8))*cos(x(7))-cos(x(9))*sin(x(7)));
xdot(6)= T/m*(cos(x(8))*cos(x(7)))-g;
xdot(7:9)= eulerrates;
xdot(10:12)=omegadot;
end
