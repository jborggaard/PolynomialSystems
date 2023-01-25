function [A,B,N] = vanderPolRing(g,a,b,c)
%-------------------------------------------------------------------------------
%vanderPolRing returns a polynomial system for a ring of van der Pol oscillators
%
%  A ring of n van der Pol oscillators with a term modeling coupled springs
%  takes the form
%
%   \ddot{x}_i = a_i (1-x_i^2) \dot{x}_i - x_i + b_i(x_{i-1}-2x_i+x_{i+1}) + u_i
%
%  where the vectors a and b are optional input vectors.  This system is studied
%  in a number of papers including:
%    Barron, Stability of a ring of coupled van der Pol oscillators with
%      non-uniform distribution of the coupling parameter, J. Appl. Res. Tech.,
%      v. 14, pp. 62--66, 2016.
%    Woafo and Enjieu, Synchronization states in a ring of mutuallycoupled
%      selfsustained electrical oscillators, Phys Rev E, v. 69, pp. 1--9, 2004.
%
%  The system is returned as a first-order system of 2*n differential equations
%  in the form
%
%    \dot{y} = A*y + B*u + N{3}*kron(kron(y,y),y)
%
%  N{3} is the only non-empty entry of the cell array N containing a matrix
%  of dimension 2*n x (2*n)^3.
%
%  Usage: [A,B,N] = vanderPolRing(n,cIndex,a,b)
%
%  by default, a = b = ones(n,1) and cIndex is a list of integers between 1 and
%  n that determine which oscillators contain a control actuator.  Example,
%
%  >> [A,B,N] = vanderPolRing(8,[1 3]);
%  >> rhs = @(t,x) A*x + N{3}*kron(kron(x,x),x);  % open loop system
%
%  Author: Jeff Borggaard
%
%  License: MIT
%
%  Part of the PolynomialSystems repository:
%          https://github.com/jborggaard/PolynomialSystems
%-------------------------------------------------------------------------------
%%
  if ( ~exist('g','var') )
    g      = 4;    % number of van der Pol oscillators
    Cidx   = [1 2];
  end

  n      = 2*g;
  m      = length(Cidx);

  if ( ~exist('a','var') )
    a  = ones(g,1);   % viscous damping parameter
  end
  if ( ~exist('b','var') )
    b  = ones(g,1);   % coupling parameter
  end

  A  = zeros(n,n);

  for i=1:g
    i1 = 2*i-1;
    i2 = 2*i;
    A(i1,i2) = 1;
    A(i2,i2) = a(i);
    A(i2,i1) = -(1+2*b(1));

    if ( i==1 )
      im1 = 2*g-1;
      im2 = 2*g;
    else
      im1 = 2*(i-1)-1;
      im2 = 2*(i-1);
    end

    if ( i==g )
      ip1 = 1;
      ip2 = 2;
    else
      ip1 = 2*(i+1)-1;
      ip2 = 2*(i+1);
    end

    A(i2,ip1) = A(i2,ip1) + b(i);
    A(i2,im1) = A(i2,im1) + b(i);
  end

  N{2} = zeros(n,n^2);
  N{3} = zeros(n,n^3);

  % set the N3 terms here...
  idx3 = @(i1,i2,i3) i1 + (i2-1)*n + (i3-1)*n^2;

  for i=1:g
    i1 = 2*i-1;
    i2 = 2*i;
    N{3}(i2,idx3(i1,i1,i2)) = N{3}(i2,idx3(i1,i1,i2)) - a(i)/3;
    N{3}(i2,idx3(i1,i2,i1)) = N{3}(i2,idx3(i1,i2,i1)) - a(i)/3;
    N{3}(i2,idx3(i2,i1,i1)) = N{3}(i2,idx3(i2,i1,i1)) - a(i)/3;
  end

  B  = zeros(n,m);

  for i=1:m
    cIndex = Cidx(i);
    B(2*cIndex,i) = 1;
  end
  
end
