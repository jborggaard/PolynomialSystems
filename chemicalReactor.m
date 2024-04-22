function [A,B,Nxx,Nxu,Q,R] = chemicalReactor()
%chemicalReactor
%  This example from
%  Hofer and Tibken, "An iterative method for the finite-time bilinear
%     quadratic control problem", JOTA, v57, pp. 411-427, 1988.
%
%  is used in the paper 
%   simulate from [.15;0] from [0,3].
%
%  setPaths
%  [k,v] = pqrBilinear(A,B,Q,R,Nxx,Nxu,zeros(2,1),3);
%  u = @(x) kronPolyEval(k,x);
%  rhs = @(t,x) A*x + B*u(x) + Nxu{1}*kron(x,u(x));
%  [T,X] = ode23(rhs,[0 3],[.15;0]);
%  plot(T,X(:,1))
%%

  A = [13/6 5/12;-50/3 -8/3];
  B = [-1/8; 0];
  Nxx{2} = zeros(2,4);
  Nxu{1} = [-1 0;0 0];

  Q = 10*eye(2);
  R = 1;

end
