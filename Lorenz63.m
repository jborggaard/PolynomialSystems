function [A,B,N] = Lorenz63(sigma,rho,beta)
%Lorenz63 returns the quadratic model for the Lorenz63 system
%   \dot{x1} = sigma*(x2-x1) + u
%   \dot{x2} = x1*(rho-x3)-x2
%   \dot{x3} = x1*x2 - beta*x3
%
%  where u is a control input.  This system is written as
%
%     \dot{x} = A*x + B*u + N*kron(x,x)
%
%  and this function returns the matrices A, B, and N (symmetrized).
%
%  Usage:  [A,B,N] = Lorenz63(sigma,rho,beta)
%
%  by default, sigma=10, beta=8/3, and rho=28 to match Lorenz' version of the
%  system from his J. Atmos. Sci. paper (v.20, pp. 130-141, 1963).
%
%  Author: Jeff Borggaard
%
%  License: MIT
%
%  Part of the PolynomialSystems repository: 
%          https://github.com/jborggaard/PolynomialSystems
%% 
  
  if (~exist('sigma','var'))
    sigma = 10;
  end

  if (~exist('rho','var'))
    rho = 28;
  end

  if (~exist('beta','var'))
    beta = 8/3;
  end
  
  A = [-sigma  sigma    0  ; ...
         rho    -1      0  ; ...
          0      0   -beta ];
  B = [1;0;0];
  N = zeros(3,9); 
  N(2,3)=-0.5; N(2,7)=-0.5;  % -x1 x3 term
  N(3,2)= 0.5; N(3,4)= 0.5;  %  x1 x2 term
  
end
