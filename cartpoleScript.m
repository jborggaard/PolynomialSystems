%  Script to run the cartpole example
%  The states are:  x(1) - cart position
%                   x(2) - cart velocity
%                   x(3) - pole angle (measured positive from vertical)
%                   x(4) - pole angular velocity
%%

addpath('/Volumes/borggaard4/Software/MyPublicSoftware/QQR')

controlDegree = 3;

g = 9.8;  parameters.Gravity = g;
m = 0.5;  parameters.MassPole = m;
M = 0.5;  parameters.MassCart = M;
L = 0.8;  parameters.Length = L;

[A,B,Nxx,Nxu] = cartpolePolynomial(m,M,L,g);

Q = ones(4);  Q(1,1) = 10;  Q(3,3) = 10;
R = 1;
parameters.Q = Q;
parameters.R = R;

%  Calculate the feedback by linearizing the model and applying LQR
[K,V] = lqr(A,B,Q,R);
[k,v] = pqrBilinear(A,B,Q,R,Nxx,Nxu,zeros(4,1),controlDegree);

%u = @(x) 0;                         % for open-loop simulation
%u = @(x) -K*x(1:4);                 % for linear feedback simulation
u = @(x) kronPolyEval(k,x(1:4));    % for cubic feedback simulation

x0 = [2.00;0.2;-0.1;0.1];           % good value function approximation
%x0 = [2.00; 0.0; -pi/3.3; 0.00];    % both feedback controls work here
%x0 = [2.00; 0.0; -pi/3  ; 0.00];    % only the cubic feedback works here

tRange = linspace(0,16,201);  % evenly spaced times makes the animation below
                              % more accurate.
[T,Y] = cartpoleIntegration(x0,tRange,u,parameters);

cartpoleAnimation(T,Y(:,1:4),parameters)

fprintf('The value function approximation is %g\n',kronPolyEval(v,x0))

