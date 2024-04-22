function [A,B,Nxx,Nxu] = acrobotPolynomial(x0,parameters)
%acrobotPolynomial  Computes a polynomial model for the acrobot using symbolics
%
%    Note: this function uses the symbolic toolbox
%
%  The acrobot model is a classical two-pendulum system that is controlled
%  through a torque at the joint between the pendula.  This mimics a gymnast on
%  the high bar, though here the joint can move 360 degrees, thus acrobot.
%
%  This function provides a cubic polynomial model expanded about a given 
%  equilibrium point x_eq, choices include:
%       [pi;0;0;0], [0;pi;0;0], [pi;pi;0;0], and [0;0;0;0].  
% 
%  The resulting model has the form:
%
%    \dot{x} = Ax + Bu + Nxx{3}*kron(kron(x,x),x) + Nxu{2}*kron(kron(x,x),u)
%
%  Usage:
%     [A,B,Nxx,Nuu] = acrobotPolynomial(x_eq,parameters)
%
%  Variables:
%     parameters - a struct with the following fields
%       m1      - mass of the leg connected to the base
%       m2      - mass of the leg that forms the second pendulum
%       l1,l2   - corresponding lengths of the legs
%       lc1,lc2 - distances along the legs to their center of mass
%       I1,I2   - corresponding rotary moments of inertia
%       g       - gravitational constant
%
%
%  Author:
%     Jeff Borggaard, Virginia Tech
%
%  Based on the model presented at underactuated.mit.edu, Chapter 3
%  by Russ Tedrake.
%
%  License: MIT
%
%  Part of the PolynomialSystems repository.
%%

  l1  = parameters.l1;
  l2  = parameters.l2;
  lc1 = parameters.lc1;
  lc2 = parameters.lc2;
  m1  = parameters.m1;
  m2  = parameters.m2;
  I1  = parameters.I1;
  I2  = parameters.I2;
  g   = parameters.g;

  %
  %  Build the right-hand-side symbolically as f(x) + g(x)u
  syms('x',[4,1])   %#ok  (we are using the variables x1, etc this generates

  M = [I1 + I2 + m2*l1^2 + 2*m2*l1*lc2*cos(x2), I2 + m2*l1*lc2*cos(x2);   ...
       I2 + m2*l1*lc2*cos(x2)                 , I2                    ];

  tau = [-m1*g*lc1*sin(x1) - m2*g*(l1*sin(x1) + lc2*sin(x1+x2)); ...
         -m2*g*lc2*sin(x1+x2)                                  ];

  C = [-2*m2*l1*lc2*sin(x2)*x4,  -m2*l1*lc2*sin(x2)*x4; ...
          m2*l1*lc2*sin(x2)*x3, 0                     ];

  f = [x3               ; ...
       x4               ; ... 
       M\(tau-C*[x3;x4])];
  g = [zeros(2,1); ...
       M\[0;1]   ];

  %
  %  Get the first-order terms for Ax + Bu
  J = jacobian(f,[x1 x2 x3 x4]);
  A = double(subs(J,{x1,x2,x3,x4},x0.'));

  B = double(subs(g,{x1,x2,x3,x4},x0.'));

  %N2 = zeros(4,16);
  x = [x1 x2 x3 x4];
  for j=1:4
    for i=1:4
      N2(:,i+4*(j-1)) = diff(J(:,i),x(j));
    end
  end
  
  N2 = N2/2; % either to remove 2 in differentiating xi^2 or account for repeats.
  
  Nxx{2} = double(subs(N2,{x1,x2,x3,x4},x0.'));  % should be zero
  
  B1 = jacobian(g,[x1 x2 x3 x4]);
  Nxu{1} = double(subs(B1,{x1,x2,x3,x4},x0.'));  % should be zero
  
  %N3 = zeros(4,64);
  for j=1:4
    for i=1:16
      N3(:,i+16*(j-1)) = diff(N2(:,i),x(j));
    end
  end
  
  N3 = N3/3;  % account for xi^3 term or repeats
  
  Nxx{3} = double(subs(N3,{x1,x2,x3,x4},x0.'));
  
  B2 = [diff(B1,x1) diff(B1,x2) diff(B1,x3) diff(B1,x4)];
  
  B2 = B2/2;  % account for xi^2 term or repeats
  
  Nxu{2} = double(subs(B2,{x1,x2,x3,x4},x0.'));

  %  For higher degree term approximations
  for j=1:4
    for i=1:64
      N4(:,i+64*(j-1)) = diff(N3(:,i),x(j));
    end
  end

  N4 = N4/4;

  Nxx{4} = double(subs(N4,{x1,x2,x3,x4},x0.'));

  B3 = [diff(B2,x1) diff(B2,x2) diff(B2,x3) diff(B2,x4)];

  B3 = B3/3;

  Nxu{3} = double(subs(B3,{x1,x2,x3,x4},x0.'));

  for j=1:4
    for i=1:256
      N5(:,i+256*(j-1)) = diff(N4(:,i),x(j));
    end
  end

  N5 = N5/5;

  Nxx{5} = double(subs(N5,{x1,x2,x3,x4},x0.'));

  B4 = [diff(B3,x1) diff(B3,x2) diff(B3,x3) diff(B3,x4)];

  B4 = B4/4;

  Nxu{4} = double(subs(B4,{x1,x2,x3,x4},x0.'));
  
end
