function cartpoleAnimation(t,x,parameters)
%cartpoleAnimation Provides an animation of the cartpole from time and states
%
%  Usage:
%       cartpoleAnimation(t,x,parameters)
%
%  Variables:
%       t  - a time vector of dimension
%       x  - a state vector of dimension
%       parameters - a struct containing fields (l1 and l2) the lengths of the
%                    two links
%
%  Adapted from a similar function, AcrobatLqr.m by Krishna Prakash Yadav 
%  on Matlab's File Exchange.
%
%  Author: Jeff Borggaard
%
%  License: MIT
%
%  Part of the PolynomialSystems repository:
%          https://github.com/jborggaard/PolynomialSystems
%%

  % cart dimensions (for visualization)
  hh = 0.05; % half height
  hw = 0.25; % half width

  % pole dimensions
  L = parameters.Length;

  % get plot range
  xmin = min(0,min(x(:,1))-hw);
  xmax = max(0,max(x(:,1))+hw);
  ymin = -hh;  % perhaps -L to account for swing up or open-loop cases.
  ymax = L;

  xmin = 1.2*xmin; xmax = 1.2*xmax;  % add a 20% border
  ymin = 1.2*ymin; ymax = 1.2*ymax;


  K = length(t);   % number of uniform timesteps in the simulation

  % initialize the cart and link positions
  x1 = x(1,1);
  y1 = 0.0;

  theta = x(1,3);
  x2 = x1 + L*sin(theta);
  y2 = y1 + L*cos(theta);

  cart = line([x1-hw x1+hw x1+hw x1-hw x1-hw],[-hh -hh hh hh -hh],...
              'color','k','LineWidth',4);
  pole = line([x1 x2],[y1 y2],...
              'color','r','LineWidth',4);

  T = t(1);
  str = strcat(num2str(T)+"/",num2str(T(end)));
  timestamp = text(1.6,ymax-0.1,str,'HorizontalAlignment','right');
  axis([xmin xmax ymin ymax])

  for i = 1:K
    theta = x(i,3);
    x1 =  x(i,1);
    y1 = 0.0;
    x2 = x1 + L*sin(theta);
    y2 = y1 + L*cos(theta);
    set(cart,'xdata',[x1-hw x1+hw x1+hw x1-hw x1-hw],'ydata',[-hh -hh hh hh -hh]);
    set(pole,'xdata',[x1 x2],'ydata',[y1 y2]);
    T = t(i);
    str = strcat(num2str(T)+"/",num2str(t(end)));
    set(timestamp,'String',str);
    pause(0.02);
    drawnow
  end
end
