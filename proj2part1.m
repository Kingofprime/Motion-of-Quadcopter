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
r0vec=[0,0,0]';
rdot0=[0,0,0]';
psi0=0;
theta0=0;
phi0=0;
p0=0;
q0=0;
r0=0;
x0=[r0vec; rdot0; psi0; theta0; phi0; p0; q0; r0];
t=6;
tspan = linspace(0,t,1000);

%% set ode45

options = odeset('AbsTol',1e-12,'RelTol',1e-12);
[t1,p1]=ode45(@(t,x)ode(t,x,g,m,l,k,b,I),tspan,x0,options);

%% Interial Frame to Body Frame
for i=1:length(tspan)
    dcm1= [1,0,0 ; 0, cos(p1(i,7)), -sin(p1(i,7)); 0 , sin(p1(i,7)), cos(p1(i,7))];
    dcm2=[cos(p1(i,8)) , 0 , sin(p1(i,8)); 0,1,0; -sin(p1(i,8)), 0 , cos(p1(i,8))];
    dcm3 = [cos(p1(i,9)) , -sin(p1(i,9)), 0; sin(p1(i,9)), cos(p1(i,9)), 0; 0,0,1];
    dcm=dcm1*dcm2*dcm3;
    bf(i,:)=dcm*[p1(i,4);p1(i,5);p1(i,6)];
end
%% Plots

% Position vs Time
figure
subplot(3,1,1)
plot(t1,p1(:,1));
ylabel('X position (m)');
xlabel('time (sec)');
title('X position in intertial frame vs time');
subplot(3,1,2)
plot(t1,p1(:,2));
ylabel('Y position (m)');
xlabel('time (sec)');
title('Y position in intertial frame vs time');
subplot(3,1,3)
plot(t1,p1(:,3));
ylabel('Z position (m)');
xlabel('time (sec)');
title('Z position in intertial frame vs time');

% Velocity vs Time
figure
subplot(3,1,1)
plot(t1,bf(:,1));
ylabel('X velocity (m\s)');
xlabel('time (sec)');
title('X velocity in body frame vs time');
subplot(3,1,2)
plot(t1,bf(:,2));
ylabel('Y velocity (m\s)');
xlabel('time (sec)');
title('Y velocity in body frame vs time');
subplot(3,1,3)
plot(t1,bf(:,3));
ylabel('Z velocity (m\s)');
xlabel('time (sec)');
title('Z velocity in body frame vs time');

% Euler angle vs Time
figure
subplot(3,1,1)
plot(t1,p1(:,7));
ylabel('phi (rad)');
xlabel('time (sec)');
title('Euler angle (phi) vs time');
subplot(3,1,2)
plot(t1,p1(:,8));
ylabel('theta (rad)');
xlabel('time (sec)');
title('Euler angle (theta) vs time');
subplot(3,1,3)
plot(t1,p1(:,9));
ylabel('psi (rad)');
xlabel('time (sec)');
title('Euler angle (psi) vs time');

% Angular velocity vs Time
figure
subplot(3,1,1)
plot(t1,p1(:,10));
ylabel('p (rad/sec)');
xlabel('time');
title('Angular Velocity (p) vs time (sec)');
subplot(3,1,2)
plot(t1,p1(:,11));
ylabel('q (rad/sec)');
xlabel('time (sec)');
title('Angular Velocity (q) vs time');
subplot(3,1,3)
plot(t1,p1(:,12));
ylabel('r (rad/sec)');
xlabel('time (sec)');
title('Angular Velocity (r) vs time');


%% Functions for ode45

function [xdot] = ode(t,x,g,m,l,k,b,I)

omhover = sqrt(m*g/(4*k));

if (t<1)
    om = [omhover + (70*sin(2*pi*t/4));omhover + (70*sin(2*pi*t/4));omhover + (70*sin(2*pi*t/4));omhover + (70*sin(2*pi*t/4))];
elseif(t>=1) && (t<=2)
    om = [omhover - (77*sin(2*pi*t/4));omhover - (77*sin(2*pi*t/4));omhover - (77*sin(2*pi*t/4));omhover - (77*sin(2*pi*t/4))];
elseif (t>=2) && (t<=3)
    om = [omhover;sqrt(omhover^2 - (70^2*sin(2*pi*(t-2)/4)));omhover;sqrt(omhover^2 + (70^2*sin(2*pi*(t-2)/4)))];
elseif (t>=3) && (t<=4)
    om = [omhover;sqrt(omhover^2 + (70^2*sin(2*pi*(t-2)/4)));omhover;sqrt(omhover^2 - (70^2*sin(2*pi*(t-2)/4)))];
elseif (t>=4) && (t<=5)
    om = [sqrt(omhover^2 - (70^2*sin(2*pi*(t-4)/4)));omhover;sqrt(omhover^2 + (70^2*sin(2*pi*(t-4)/4)));omhover];
elseif (t>=5) && (t<=6)
    om = [sqrt(omhover^2 + (70^2*sin(2*pi*(t-4)/4)));omhover;sqrt(omhover^2 - (70^2*sin(2*pi*(t-4)/4)));omhover];
end

L=k*l*(-om(2)^2+om(4)^2);
M=k*l*(-om(1)^2+om(3)^2);
N=b*(om(1)^2-om(2)^2+om(3)^2-om(4)^2);
tau = [L;M;N];

omega = x(10:12);
omegadot = I\(cross(-omega,I*omega)+tau);

A= [1,0, -sin(x(8)) ; 0, cos(x(7)), cos(x(8))*sin(x(7)); 0, -sin(x(7)), cos(x(8))*cos(x(7))];

eulerrates=A\x(10:12);
T = k*sum(om.^2);

% output
xdot=zeros(size(x));
xdot(1:3)=x(4:6);
xdot(4)= T/m*(cos(x(9))*sin(x(8))*cos(x(7))+sin(x(9))*sin(x(7)));
xdot(5)= T/m*(sin(x(9))*sin(x(8))*cos(x(7))-cos(x(9))*sin(x(7)));
xdot(6)= T/m*(cos(x(8))*cos(x(7)))-g;
xdot(7:9)= eulerrates;
xdot(10:12)=omegadot;
end

