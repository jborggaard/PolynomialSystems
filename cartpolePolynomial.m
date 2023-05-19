function [A,B,Nxx,Nxu] = cartpolePolynomial(m,M,L,g)
%cartpoleP  Produces a polynomial model of the frictionless cart-pole system
%  The cart-pole model seeks to control a cart that holds an inverted 
%  pendulum.  The cart has mass M, the pole has mass m, and length L,
%  while g is the gravitational constant.  The full-order model has the
%  form
%      \dot{x} = [ x(2);
%                  r*l*x(4)^2*sin(x(3)) -
%                  r*l*cos(x(3))/(l*(4/3)-r*cos(x(3))^2)*(g*sin(x(3) -
%                  l*r*x(4)^2*sin(2*x(3))/2;
%                  x(4);
%                  (g*sin(x(3))-l*r*x(4)^2*sin(2*x(3))/2)/(l*(4/3-r*cos^2(x(3)))
%                ] + B(x)*u
%   
%  This function returns a polynomial approximation to this system of the
%  form:
%      \dot{x} = A*x + B*u + Nxx{3}*kron(x,kron(x,x)) + 
%                Nxx{5}*kron(x,kron(x,kron(x,kron(x,x)))) + 
%                Nxu{2}*kron(x,x)*u
%
%  for the purpose of developing a polynomial feedback law.
%
%  Usage:
%     [A,B,Nxx,Nuu] = cartpolePolynomial(m,M,L,g)
%
%  Variables:
%     m - mass of the pole
%     M - mass of the cart
%     L - length of the pole
%     g - gravitational constant
%
%
%  Author: Jeff Borggaard
%
%  License: MIT
%
%  Part of the PolynomialSystems repository.
%%
  r = m/(m+M);

  A = [0 1 0 0; ...
       0 0 g*r/(r-4/3) 0; ...
       0 0 0 1; ...
       0 0 g/(L*(4/3-r)) 0];

  B = [0; 4*r/(3*m*(4/3-r)); 0; -r/(m*L*(4/3-r))];

  Nxx{2} = zeros(4,16);

  Nxx{3} = zeros(4,64);
  Nxx{3}(2,43) = (1/6)*(16/3*g*r+2*g*r^2)/(4/3-r)^2;
  Nxx{3}(2,63) = (1/6)*(8/3*r*L)/(4/3-r);
  Nxx{3}(2,60) = Nxx{3}(2,63);
  Nxx{3}(2,48) = Nxx{3}(2,63);

  Nxx{3}(4,43) = -(1/6)*(4/3*g+5*g*r)/(L*(4/3-r)^2);
  Nxx{3}(4,63) = -(1/3)*r/(4/3-r);
  Nxx{3}(4,60) = Nxx{3}(4,63);
  Nxx{3}(4,48) = Nxx{3}(4,63);

  Nxx{4} = zeros(4,256);

  Nxx{5} = zeros(4,1024);
  % x3^5 term
  Nxx{5}(2, 683) = (1/120)*(-(256/9)*g*r - 352/3*g*r^2 - 16*g*r^3)/(4/3-r)^3;
  % x3^3 x4^2 terms
  Nxx{5}(2,1003) = (1/120)*(-(32/9)*r*L - (40/3)*L*r^2)/(4/3-r)^2;
  Nxx{5}(2, 955) = Nxx{5}(2,1003);
  Nxx{5}(2, 763) = Nxx{5}(2,1003);
  Nxx{5}(2, 943) = Nxx{5}(2,1003);
  Nxx{5}(2, 751) = Nxx{5}(2,1003);
  Nxx{5}(2, 703) = Nxx{5}(2,1003);
  Nxx{5}(2, 940) = Nxx{5}(2,1003);
  Nxx{5}(2, 748) = Nxx{5}(2,1003);
  Nxx{5}(2, 700) = Nxx{5}(2,1003);
  Nxx{5}(2, 688) = Nxx{5}(2,1003);

  Nxx{5}(4, 683) = (1/120)*((16/9)*g + (232/3)*g*r + 61*g*r^2)/(L*(4/3-r)^3);

  Nxx{5}(4,1003) = (1/120)*(4*r^2+(32/3)*r)/(4/3-r)^2;
  Nxx{5}(4, 955) = Nxx{5}(4,1003);
  Nxx{5}(4, 763) = Nxx{5}(4,1003);
  Nxx{5}(4, 943) = Nxx{5}(4,1003);
  Nxx{5}(4, 751) = Nxx{5}(4,1003);
  Nxx{5}(4, 703) = Nxx{5}(4,1003);
  Nxx{5}(4, 940) = Nxx{5}(4,1003);
  Nxx{5}(4, 748) = Nxx{5}(4,1003);
  Nxx{5}(4, 700) = Nxx{5}(4,1003);
  Nxx{5}(4, 688) = Nxx{5}(4,1003);

  Nxu{1} = zeros(4,4);
  
  Nxu{2} = zeros(4,16);
  Nxu{2}(2,11) = -r^2/(m*(4/3-r)) -r^3/(m*(4/3-r)^2);
  Nxu{2}(4,11) = (1/2)*r/(m*L*(4/3-r)) + r^2/(m*L*(4/3-r)^2);
  
end