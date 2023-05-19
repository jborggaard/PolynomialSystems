function acrobotAnimation(t,x,parameters)
  %acrobotAnimation Provides an animation of the acrobot from time and states
  %
  %  Usage:
  %       acrobotAnimation(t,x,parameters)
  %
  %  Variables:
  %       t  - a time vector of dimension
  %       x  - a state vector of dimension
  %       parameters - a struct containing fields (l1 and l2) the lengths of the
  %                    two links
  %
  %  Adapted from AcrobatLqr.m 
  %
  %  Part of the PolynomialSystems repository.
  %%

  l1 = parameters.l1;
  l2 = parameters.l2;
 
  K = length(t);   % number of uniform timesteps in the simulation

  % initialize the end positions of the two links
  theta1 = x(1,1);
  x1 =  l1*sin(theta1);
  y1 = -l1*cos(theta1);

  theta2 = x(1,2);
  x2 = x1 + l2*sin(theta1+theta2);
  y2 = y1 - l2*cos(theta1+theta2);

  link1 = line([ 0 x1],[ 0 y1],'color','k','LineWidth',4);
  link2 = line([x1 x2],[y1 y2],'color','r','LineWidth',4);

  T = t(1);
  str = strcat(num2str(T)+"/",num2str(T(end)));
  timestamp = text(1.6,-0.1,str,'HorizontalAlignment','right');
  axis([-3 3 -1 3])

  for i = 1:K
    theta1 = x(i,1);
    theta2 = x(i,2);
    x1 =  l1*sin(theta1);
    y1 = -l1*cos(theta1);
    x2 = x1 + l2*sin(theta1+theta2);
    y2 = y1 - l2*cos(theta1+theta2);
    set(link1,'xdata',[ 0 x1],'ydata',[ 0 y1]);
    set(link2,'xdata',[x1 x2],'ydata',[y1 y2]);
    T = t(i);
    str = strcat(num2str(T)+"/",num2str(t(end)));
    set(timestamp,'String',str);
    pause(0.02);
    drawnow
  end
end