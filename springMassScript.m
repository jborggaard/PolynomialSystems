%  This script evaluates the performance of polynomial feedback control for
%  a spring-mass system with a hardening/softening spring model.
%
%    m \ddot{x} = -c \dot{x} - k x - beta x^3 + u(t)       (eq. 1)
% 
%    m    - mass
%    c    - damping
%    k    - stiffness (Hooke's Law term)
%    beta - nonlinear stiffness coefficient  
%           (beta>0 is hardening, beta<0 is softening)
%
%  Note:  To turn this example into a simple pendulum, 
%
%           \ddot{theta} = -g/L sin(\theta) + u(t)
%
%         set c=0, k= g/L, beta=-g/(6L) where L is the length of the pendulum
%         and g is the gravitational constant (the mass cancels so, set m=1 to%         have the form of the spring model above)
%
%  Part of the PolynomialSystems repository.
%%

addpath('../QQR')
addpath( '../KroneckerTools/src' )

%  Set physical parameters
m    = 1;
c    = 0;
ks   = 1;
beta = -1/6;

% %  Here are parameters for a cubic model of the simple pendulum
% m    = 1;
% c    = 0;
% L    = 1;
% g    = 9.81;
% k    = g/L;
% beta =-g/(6*L);

%  Set control parameters
Q    = [1 0;0 0];  % we are only explicitly controlling the displacement
R    = .1 ;        % making this smaller gives the control more "authority"

%  Setup inputs to QQR (writing as a first-order system with polynomial RHS)
A    = [0 1;-ks/m -c/m];
B    = [0;1];
N{2} = zeros(2,4);
N{3} = zeros(2,8);
N{3}(2,1) = -beta/m;   % the -beta/m x(1)^3 term in (eq. 1)

% compute polynomial approximations to the feedback (k) and value function (v)
degree = 5;
[k,v] = pqr(A,B,Q,R,N,degree);

% set the initial condition to test control effectiveness
x0 = [2;0];    
tspan = [0 10];

% the open-loop model without any control force
rhs_open = @(t,x) A*x + N{3}*kron(kron(x,x),x);

% try different degrees of feedback control
u1 = @(x) kronPolyEval(k,x,1);
u3 = @(x) kronPolyEval(k,x,3);
u5 = @(x) kronPolyEval(k,x,5);
rhs_1 = @(t,x) A*x + N{3}*kron(kron(x,x),x) + B*u1(x);
rhs_3 = @(t,x) A*x + N{3}*kron(kron(x,x),x) + B*u3(x);
rhs_5 = @(t,x) A*x + N{3}*kron(kron(x,x),x) + B*u5(x);

% simulate the open-loop
[T0,X0] = ode45(rhs_open,tspan,x0);

% simulate the closed-loop cases
[T1,X1] = ode45(rhs_1,tspan,x0);
[T3,X3] = ode45(rhs_3,tspan,x0);
[T5,X5] = ode45(rhs_5,tspan,x0);

plot(T0,X0(:,1),'r',...
     T1,X1(:,1),'g*-',...
     T3,X3(:,1),'b+-',...
     T5,X5(:,1),'k-')
xlabel('time'); ylabel('displacement'); 
title('Nonlinear Spring-Mass with Control')
legend('no control',...
       'degree 1 feedback',...
       'degree 3 feedback',...
       'degree 5 feedback')


%  Now let's repeat the control experiment, but also estimate the value 
%  function:  int_0^T x'*Q*x + u'*R*u dt, we do this by augmenting the
%  closed-loop ODEs with another equation:  c_dot = x'*Q*x + u'*R*u that is
%  integrated from c(0)=0 to c(T).  We use a "trim function" to pull out
%  the real state from x since the last entry of x is now c.  Since the value %  function is the integral with T=\infty, we increase the final time in this
%  approximation.
tspan = [0 20];
trim = @(x) x(1:end-1);
rhs_1E = @(t,x) [A*trim(x) + N{3}*kron(kron(trim(x),trim(x)),trim(x)) + B*u1(trim(x)); ...
                trim(x)'*Q*trim(x) + u1(trim(x))'*R*u1(trim(x))];
[~,X1] = ode45(rhs_1E,tspan,[x0;0]);
fprintf('For degree 1 feedback, degree 2 value function\n')
fprintf('Approximate value function: %g,  Integrated running cost: %g\n\n',...
        kronPolyEval(v,x0,2),    X1(end,end))

rhs_3E = @(t,x) [A*trim(x) + N{3}*kron(kron(trim(x),trim(x)),trim(x)) + B*u3(trim(x)); ...
                 trim(x)'*Q*trim(x) + u3(trim(x))'*R*u3(trim(x))];
[~,X3] = ode45(rhs_3E,tspan,[x0;0]);
fprintf('For degree 3 feedback, degree 4 value function\n')
fprintf('Approximate value function: %g,  Integrated running cost: %g\n\n',...
        kronPolyEval(v,x0,4),    X3(end,end))

rhs_5E = @(t,x) [A*trim(x) + N{3}*kron(kron(trim(x),trim(x)),trim(x)) + B*u5(trim(x)); ...
                 trim(x)'*Q*trim(x) + u5(trim(x))'*R*u5(trim(x))];
[~,X5] = ode45(rhs_5E,tspan,[x0;0]);
fprintf('For degree 5 feedback, degree 6 value function\n')
fprintf('Approximate value function: %g,  Integrated running cost: %g\n\n',...
        kronPolyEval(v,x0,6),    X5(end,end))
