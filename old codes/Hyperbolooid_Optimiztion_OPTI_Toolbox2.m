clc;
clear;
close all;

Vmin = 1.820219934;
Smax = 4.523893422;
hend = 3*Vmin/Smax;

i = 0;
XX = [];
for h = 0.9:0.1:hend
    i = i + 1;
    hh(i) = h;
% Objective
% x1 = a, x2 = c
fun = @(x) (576./((pi.*h./3)*(3*x(2)^2 + h^2).*(x(1)./x(2)).^2)).*(pi/2).*(x(1)./x(2)).^4.*(h.^5./5 + (2.*x(2).^2.*h.^3)./3 + x(2).^4.*h);     
% fun = @(x) (576./((pi./3).*(x(2)-h).^2.*(2.*x(2)+h).*(x(1)./x(2)).^2)).*(pi./2).*(x(1)./x(2)).^4.*((h.^5-x(2).^5)./5-2.*x(2).^2.*(h.^3-x(2).^3)./3+x(2).^4.*(h-x(2)));

% Bounds (lb <= x <= ub)
A = (sqrt(Vmin*h*(-Smax*h + 3*Vmin)^2/((Smax*h^3 + Vmin*h^2 + h*(-Smax*h + 3*Vmin)^2)*h*pi)));
C = (sqrt(h^2*(-Smax*h + 3*Vmin)/(3*(Smax*h + Vmin))));
% A = sqrt(6.*(Smax.*h - 3.*Vmin + sqrt(9.*Smax.^2*h.^2 - 30.*Smax.*Vmin.*h + 9.*Vmin.^2)).^2.*Vmin.*Smax./(pi*(3.*Smax.*h - 3.*Vmin + sqrt(9.*Smax.^2.*h.^2 - 30.*Smax.*Vmin.*h + 9.*Vmin.^2)).*(3.*Smax.*h + 3.*Vmin - sqrt(9.*Smax.^2.*h.^2 - 30.*Smax.*Vmin.*h + 9.*Vmin.^2)).^2));
% C = (Smax.*h - 3.*Vmin + sqrt(9.*Smax.^2*h.^2 - 30.*Smax.*Vmin.*h + 9.*Vmin.^2))./(4.*Smax);
lb = [0, 0];
ub = [A, C]; 

% Nonlinear Constraints (cl <= nlcon(x) <= cu)
% nlcon = @(x) pi*(x(1)/x(2))^2*(h^2-x(2)^2);
% cl = 1;
% cu = 4.5239; 

% Initial Guess
x0 = [1 1]';

% Options
opts = optiset('solver','ipopt','display','iter');

% Build OPTI Problem
Opt = opti('fun',fun,'bounds',lb,ub,'x0',x0,'options',opts);

% Solve NLP
[x,fval,exitflag,info] = solve(Opt);

figure
plot(Opt)

% Checking Volume and Area
% Vol(i) = (1/3)*pi*h*(x(1)/x(2))^2*(3*x(2)^2+h^2);
% Are(i) = pi*(x(1)/x(2))^2*(h^2-x(2)^2);

Vol(i) = (pi./3).*(x(2)-h).^2.*(2.*x(2)+h).*(x(1)./x(2)).^2;
Are(i) = pi*(x(1)./x(2)).^2.*(h.^2-x(2).^2);

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

figure(6);
plot(XX(1,:),XX(2,:))
xlabel('X1'); ylabel('X2')

