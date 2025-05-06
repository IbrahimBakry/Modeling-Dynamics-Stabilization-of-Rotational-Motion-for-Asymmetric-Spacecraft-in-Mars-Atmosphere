% Objective
% x1 = h, x2 = a, x3 = c
fun = @(x) (3*576/(pi*1.02*(3*x(3)^2 + 1.02^2)*(x(2)/x(3))^2))*(pi/2)*(x(2)/x(3))^4*(x(1)^5/5 + (2/3)*(x(3)^2*x(1)^3) + x(3)^4*x(1));

% Nonlinear Constraints (cl <= nlcon(x) <= cu)
nlcon = @(x) [ (1/3)*pi*x(1)*(x(2)/x(3))^2*(3*x(3)^2+x(1)^2);
               pi*(x(2)/x(3))^2*(x(1)^2-x(3)^2) ];
cl = [1.8202; 
      2.2620];
cu = [Inf; 
      4.5239];       

  
% Initial Guess
x0 = [1.02 1 1]';

% Bounds (lb <= x <= ub)
lb = [1, 1, 1];
ub = [3.0, 3.0, 3.0];  

% Options
opts = optiset('solver','ipopt','display','iter');

% Build OPTI Problem
Opt = opti('fun',fun,'nl',nlcon,cl,cu,'bounds',lb,ub,'x0',x0,'options',opts)

% Solve NLP
[x,fval,exitflag,info] = solve(Opt)

