function [T,X] = acrobotIntegration(x0,tRange,u,parameters)
%acrobotIntegration Numerically integrates the ODEs for the acrobot system.
%  These equations are consistent with the model appearing on the course webpage
%  of Russ Tedrake (Underactuated Robotics course notes).
%
%  Usage:
%      [T,X] = acrobotIntegration(x0,tRange,u,parameters)
%%
  I1 = parameters.I1;
  I2 = parameters.I2;
  l1 = parameters.l1;
  l2 = parameters.l2;
  m1 = parameters.m1;
  m2 = parameters.m2;
  g  = parameters.g;
  Q  = parameters.Q;
  R  = parameters.R;

  xe = [pi;0;0;0];
  
  M = @(x)[I1 + I2 + m2*l1^2 + m2*l1*l2*cos(x(2)), I2 + m2*l1*l2/2*cos(x(2));   ...
           I2 + m2*l1*l2/2*cos(x(2))             , I2                       ];

  tau = @(x) [-m1*g*l1/2*sin(x(1)) - m2*g*(l1*sin(x(1)) + l2/2*sin(x(1)+x(2))); ...
              -m2*g*l2/2*sin(x(1)+x(2))                                       ];

  C = @(x) [-m2*l1*l2*sin(x(2))*x(4),  -m2*l1*l2/2*sin(x(2))*x(4); ...
            m2*l1*l2/2*sin(x(2))*x(3), 0                         ];

  f = @(x) [x(3);
            x(4);
            M(x)\(tau(x)-C(x)*[x(3);x(4)])];
  b = @(x) [zeros(2,1);
            M(x)\[0;1]];

  rhs = @(t,x) [f(x(1:4)); ...
                (x(1:4)-xe).'*Q*(x(1:4)-xe) + u(x(1:4)).'*R*u(x(1:4))] + ...
               [b(x(1:4)); 0]*u(x(1:4));
  [T,X] = ode23s(rhs,tRange,[x0;0]);

  fprintf('The integrated running cost is %g\n',X(end,end))

end
