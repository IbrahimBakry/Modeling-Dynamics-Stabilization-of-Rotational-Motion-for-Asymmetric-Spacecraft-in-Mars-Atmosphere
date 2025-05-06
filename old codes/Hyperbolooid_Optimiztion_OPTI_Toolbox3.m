clc;
clear;
close all;

i = 0;
XX = [];
for h = 0.1:0.1:3
    i = i + 1;
    hh(i) = h;
% Objective
% x1 = a, x2 = c
fun = @(x) (3*576/(pi*1.02*(3*x(2)^2 + 1.02^2)*(x(1)/x(2))^2))*(pi/2)*(x(1)/x(2))^4*(h^5/5 + (2/3)*(x(2)^2*h^3) + x(2)^4*h);     

% Nonlinear Constraints (cl <= nlcon(x) <= cu)
nlcon = @(x) [ (1/3)*pi*h*(x(1)/x(2))^2*(3*x(2)^2+h^2);
               pi*(x(1)/x(2))^2*(h^2-x(2)^2) ];
cl = [1.8202; 
      2.2620];
cu = [Inf; 
      4.5239]; 

% Initial Guess
x0 = [1 1]';

% Options
opts = optiset('solver','ipopt','display','iter');

% Build OPTI Problem
Opt = opti('fun',fun,'nl',nlcon,cl,cu,'x0',x0,'options',opts);

% Solve NLP
[x,fval,exitflag,info] = solve(Opt);

% Checking Volume and Area
Vol(i) = (1/3)*pi*h*(x(1)/x(2))^2*(3*x(2)^2+h^2);
Are(i) = pi*(x(1)/x(2))^2*(h^2-x(2)^2);

% Checking the Iz
Iz(i) = fun(x);

XX(:,i) = x;
end

figure(1); 
plot(hh,XX(1,:)); 
ylabel('A'); xlabel('h')
figure(2); 
plot(hh,XX(2,:)); 
ylabel('C'); xlabel('h')

figure(3); 
plot(hh,Vol,hh,1.8202*ones(1,length(hh)))
xlabel('h'); ylabel('Volume')
figure(4); 
plot(hh,Are,hh,4.52389*ones(1,length(hh)))
xlabel('h'); ylabel('Area')

figure(5);
plot(hh,Iz)
xlabel('h'); ylabel('Iz')
