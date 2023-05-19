%  Script to run the acrobot example
addpath('/Volumes/borggaard/Software/MyPublicSoftware/QQR')

g = 9.81;  parameters.g  = g;
m1 = 1.0;  parameters.m1 = m1;
m2 = 1.0;  parameters.m2 = m2;
l1 = 1.0;  parameters.l1 = l1;   parameters.lc1 = l1/2;
l2 = 1.0;  parameters.l2 = l2;   parameters.lc2 = l2/2;
I1 = 1.0;  parameters.I1 = I1;
I2 = 1.0;  parameters.I2 = I2;

% the unknowns are 
%    x(1) - the angle of the first link to downward position
%    x(2) - the angle of the second link to the first link
%    x(3) - the time derivative of x(1)
%    x(4) - the time derivative of x(2)

% the following x0/x_eq pairs are INSIDE the LQR controlled region of attraction.
%x0 = [pi+0.05; 0.001; 0.00; 0.00];  % IF x_eq = [pi;0;0;0]
x0 = [pi/2; 0; -pi/4; pi/4];         % IF x_eq = [0;pi;0;0]

% the following x0 is OUTSIDE the LQR controlled region of attraction for xe.
%x0 = [pi+0.1; 0.0  ; 0.10; 0.00];


%xe = [pi;0;0;0];  % the acrobotPolynomial function is linearized about this
                  % (unstable) equilibrium point.
xe = [0;pi;0;0];  % the acrobotPolynomial function is linearized about this
                  % (unstable) equilibrium point.

[A,B,Nxx,Nxu] = acrobotPolynomial(xe,parameters);

Q = ones(4);  Q(1,1) = 10;  Q(3,3) = 10;
R = 0.01;
parameters.Q = Q;
parameters.R = R;

%  Calculate the feedback by linearizing the model and applying LQR
[K,V] = lqr(A,B,Q,R);
[k,v] = pqrBilinear(A,B,Q,R,Nxx,Nxu,zeros(4,1),3);

% the control is the torque applied between links 1 and 2.
%u = @(x) 0;                % for open-loop simulation
%u = @(x) -K*(x(1:4)-xe);   % for linear feedback simulation
u = @(x) kronPolyEval(k,x(1:4)-xe);

%% Simulate the closed-loop system
tRange = linspace(0,5,501);  % evenly spaced times makes the animation below
                             % more accurate.
[T,X] = acrobotIntegration(x0,tRange,u,parameters);

%% Animate the acrobot
acrobotAnimation(T,X(:,1:4),parameters)
