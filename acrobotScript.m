%  Script to run the acrobot example
%  The main parameter to change is controlDegree.  Try
%  controlDegree = 1, 3, 5 to see the benefits on nonlinear feedback.

setPaths  % find QQR and KroneckerTools

%  Control options
%    controlDegree = 0;  % open loop
%    controlDegree = 1;  % LQR based on the linearized model
%    controlDegree = 3,5 % PQR bilinear based on polynomial approximation
controlDegree = 5;

% the unknowns are 
%    x(1) - the angle of the first link to downward position
%    x(2) - the angle of the second link to the first link
%    x(3) - the time derivative of x(1)
%    x(4) - the time derivative of x(2)

%  The polynomial model is built from linearizing about the equilibrium pt: xe
xe = [pi;0;0;0];  % the acrobotPolynomial function is linearized about this
                  % (unstable) equilibrium point. (acrobot "handstand")
%xe = [0;pi;0;0];  % the acrobotPolynomial function is linearized about this
                  % (unstable) equilibrium point. (acrobot "jackknife")

%  Set the initial conditions for the open- or closed-loop simulation
% the following x0/x_eq pairs are INSIDE the LQR controlled region of attraction.
%x0 = [pi+0.1; 0.001; 0.00; 0.00];    % IF x_eq = [pi;0;0;0]
x0 = [pi+0.1; 0.001;-0.02; 0.00];    % IF x_eq = [pi;0;0;0]
%x0 = [pi/2; 0; -pi/4; pi/4];         % IF x_eq = [0;pi;0;0]

% the following x0/x_eq pairs are INSIDE degree 1 and 5 feedback, but not degree
% 3
% x0 = [pi+0.1; 0.001; 0.03; 0.00];     % IF x_eq = [pi;0;0;0]

% the following x0 is OUTSIDE the LQR controlled region of attraction for xe.
%x0 = [pi+0.1; 0; 0; pi/6];     % IF x_eq = [pi;0;0;0]



%%  Set physical parameters
g = 9.81;  parameters.g  = g;
m1 = 1.0;  parameters.m1 = m1;
m2 = 1.0;  parameters.m2 = m2;
l1 = 1.0;  parameters.l1 = l1;   parameters.lc1 = l1/2;
l2 = 1.0;  parameters.l2 = l2;   parameters.lc2 = l2/2;
I1 = 1.0;  parameters.I1 = I1;
I2 = 1.0;  parameters.I2 = I2;

%  Set control parameters
Q = ones(4);  Q(1,1) = 10;  Q(3,3) = 10;
R = 0.005;
parameters.Q = Q;
parameters.R = R;

%  Build the polynomial approximation to the acrobot
[A,B,Nxx,Nxu] = acrobotPolynomial(xe,parameters);

if ( controlDegree==0 )
   u = @(x) 0;                % for open-loop simulation

elseif ( controlDegree==1 )
   % the control is the torque applied between links 1 and 2.
   %  Calculate the feedback by linearizing the model and applying LQR
   [K,V] = lqr(A,B,Q,R);
   u = @(x) -K*(x(1:4)-xe);   % for linear feedback simulation

  vfApprox = (x0-xe).'*V*(x0-xe);
  fprintf('The degree %d approx to v(x0) is %g\n',2,vfApprox)

elseif ( controlDegree>1 )
  %  Calculate the polynomial feedback 
  [k,v] = pqrBilinear(A,B,Q,R,Nxx,Nxu,zeros(4,1),controlDegree);
  u = @(x) kronPolyEval(k,x(1:4)-xe);

  vfApprox = kronPolyEval(v,x0-xe,controlDegree+1);
  fprintf('The degree %d approx to v(x0) is %g\n',controlDegree+1,vfApprox)
end


%% Simulate the open- or closed-loop system
tRange = linspace(0,10,501);  % evenly spaced times makes the animation below
                             % more accurate.
[T,X] = acrobotIntegration(x0,tRange,u,parameters);


%% Animate the acrobot
acrobotAnimation(T,X(:,1:4),parameters)
