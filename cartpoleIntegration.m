function [T,Y] = cartpoleIntegration(x0,tRange,u,parameters)
%cartPoleIntegration Numerically integrates the ODEs for the cart-pole system.
%  These equations are consistent with the model appearing in the document
%    Florian, Correct equations for the dynamics of the cart-pole system, 2007.
%
%  For the frictionless case, this matches the classical paper of
%    Barto, Sutton, and Anderson, 1983.
%%
  m = parameters.MassPole;
  M = parameters.MassCart;
  L = parameters.Length;
  g = parameters.Gravity;
  Q = parameters.Q;
  R = parameters.R;

  r = m/(m+M);
  
  rhs = @(t,x) [x(2); ...
                r*L*x(4)^2*sin(x(3)) - r*cos(x(3))*(g*sin(x(3))-L*r*x(4)^2*sin(2*x(3))/2)/(4/3-r*cos(x(3))^2); ... 
                x(4); ...
                (g*sin(x(3))-L*r*x(4)^2*sin(2*x(3))/2)/(L*(4/3-r*cos(x(3))^2)); ...
                x(1:4).'*Q*x(1:4) + u(x(1:4))*R*u(x(1:4))] + ...
                [0; r^2*cos(x(3))^2/(m*(4/3-r*cos(x(3))^2)) + r/m; ...
                 0; -r*cos(x(3))/(m*L*(4/3-r*cos(x(3))^2)); ...
                 0]*u(x(1:4));
  
  [T,Y] = ode23s(rhs,tRange,[x0;0]);

  fprintf('The integrated running cost is %g\n',Y(end,end))

end
