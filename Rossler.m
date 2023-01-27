function [A,B,N] = Rossler(a,c)
%Rossler returns the quadratic model for the Rossler system
%   \dot{x1} =       -x2-x3
%   \dot{x2} = x1 + a*x2
%   \dot{x3} = x3*(x1-c) + u
%
%  where u is a control input.  This system is written as
%
%     \dot{x} = A*x + B*u + N*kron(x,x)
%
%  and this function returns the matrices A, B, and N (symmetrized).
%
%  Usage:  [A,B,N] = Rossler(a,c)
%
%  by default, a=0.1 and c=14 (setting u=0.1), which are commonly used choices.
%  The original parameters by Rossler were (0.2,5.7,u=0.2) (Wikipedia).
%
%  Author: Jeff Borggaard
%
%  License: MIT
%
%  Part of the PolynomialSystems repository: 
%          https://github.com/jborggaard/PolynomialSystems
%% 
  
  if (~exist('a','var'))
    a = 0.1;
  end

  if (~exist('c','var'))
    c = 5.7;
  end
  
  A = [  0    -1   -1  ; ...
         1     a    0  ; ...
         0     0   -c ];
  B = [0;0;1];
  N = zeros(3,9); 
  N(3,3)= 0.5; N(3,7)= 0.5;  %  x1 x3 term
  
end
