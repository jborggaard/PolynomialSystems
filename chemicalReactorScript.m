% for x0 = [.05;0], we are within the radius of convergence of the polynomial
% approximations to the optimal feedback and the value function.
%
% for x0 = [.1;0], the polynomial feedback laws outperform the linear feedback
% and the value function approximations get better for even polynomial
% approximations (feedbackDegree = odd).
%
% for x0 = [.15;0], the performance of linear feedback seems to be superior out
% to feedbackDegree=5.
%%
  setPaths
  x0 = [.1;0];
  feedbackDegree = 5;
 

  [A,B,Nxx,Nxu,Q,R] = chemicalReactor();
  Nuu = zeros(2,1);
  [k,v] = pqrBilinear(A,B,Q,R,Nxx,Nxu,Nuu,feedbackDegree);

  u = @(x) kronPolyEval(k,x);
  rhs = @(t,x) [A*x(1:end-1) + B*u(x(1:end-1)) + Nxu{1}*kron(x(1:end-1),u(x(1:end-1)));...
                x(1:end-1).'*Q*x(1:end-1) + u(x(1:end-1)).'*R*u(x(1:end-1))];
  [T,X] = ode23(rhs,[0 25],[x0;0]);
  
  subplot(3,1,1)
  plot(T,X(:,1))
  title('Component 1')
  subplot(3,1,2)
  plot(T,X(:,2))
  title('Component 2')
  fprintf('For feedback degree %d, the LQR cost is %g\n',feedbackDegree,X(end,end));

  nT = length(T);
  control = zeros(nT,1);
  for k=1:nT
    control(k) = u(X(k,1:end-1).');
  end
  subplot(3,1,3)
  plot(T,control)
  title('Control Input')

  fprintf('The predicted LQR cost is %g\n',kronPolyEval(v,x0,feedbackDegree+1))