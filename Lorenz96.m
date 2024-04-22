function [A,B,N] = Lorenz96(n)
%Lorenz96 returns the quadratic model for the Lorenz96 system of size n (>3).
%
%  \dot{x_1} = (x_2-x_{n-1})x_n     - x_1 + u
%  \dot{x_2} = (x_3-x_{n  })x_1     - x_2 + u
%  \dot{x_3} = (x_4-x_1    )x_2     - x_3 + u 
%    \vdots  =   \vdots
%  \dot{x_n} = (x_1-x_{n-2})x_{n-1} - x_n + u
%
%  where u is a control input.  This system is written as
%
%     \dot{x} = A*x + B*u + N*kron(x,x)
%
%  and this function returns the matrices A, B, and N (symmetrized).
%
%  Usage:  [A,B,N] = Lorenz96(n)
%
%  by default, n=4 to match Lorenz' version of the
%  system from his J. Atmos. Sci. paper (v.20, pp. 130-141, 1963).
%
%  Author: Jeff Borggaard
%
%  License: MIT
%
%  Part of the PolynomialSystems repository: 
%          https://github.com/jborggaard/PolynomialSystems
%% 
  
  %
  % Check inputs
  if (~exist('n','var'))
    n = 4;
  end

  attributes = {'>',3};
  validateattributes(n,{'numeric'},attributes)
  
  %
  % Define useful functions
  idx1 = @(i1,i2) 1+mod(i1-1,n)+mod(i2-1,n)*n;
  idx2 = @(i1,i2) idx1(i2,i1);

  %
  % Specify the system
  A = -eye(n);
  B = ones(n,1);

  N = zeros(n,n^2);
  for i=1:n
    N(i,idx1(i+1,i-1)) =  0.5;  N(i,idx2(i+1,i-1)) =  0.5;
    N(i,idx1(i-2,i-1)) = -0.5;  N(i,idx2(i-2,i-1)) = -0.5;
  end

end
