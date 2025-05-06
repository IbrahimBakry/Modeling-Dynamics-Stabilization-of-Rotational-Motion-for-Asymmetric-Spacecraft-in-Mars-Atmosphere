% Objective, a=x1, c=x2
clc;
clear;
close all;

% Mars Polar Lander: hmin=1.5, Rmax=1.25; 
% Mars3: hmin = 1.8, Rmax = 1.6;
% Insight: hmin = 1.8, Rmax = 1.3;
hmin = 1.8; % The given minimum height of KA/
Rmax = 1.6;% The final posible height of KA/
Vmin = (pi*hmin/3)*(Rmax^2);
Smax = pi*Rmax^2;

hstart= hmin;    % The given minimum height of KA/
hend = hmin+0.5; % The final posible height of KA/
t0 = 0;
tf = 300;
tstep= 20; 

k = 0;
for t = t0:tstep:tf
    k = k + 1;
[h,h1] = hlaw(hstart,hend,t,t0,tf);

m = 10; % Mass pf the deformable hyperboloid
fun = @(x) (m./(pi.*h.*(3.*x(2).^2+h.^2).*((x(1)./x(2)).^2)./3)).*(pi/2).*(x(1)./x(2)).^4.*(h.^5./5 + (2.*x(2).^2.*h.^3)./3 + x(2).^4.*h);

% Nonlinear Constraints (cl <= nlcon(x) <= cu)
nlcon = @(x) [ (1/3).*pi.*h.*(x(1)./x(2)).^2.*(3*x(2).^2+h.^2);
%                pi*(x(1)./x(2)).^2.*(h^2 + x(2).^2)]; 
                 pi*(x(1)./x(2)).^2.*(x(2).^2)]; % Area@z=h=0

cl = [Vmin;0];
cu = [Inf;Smax];

% Bounds (lb <= x <= ub)
lb = 0.1*ones(2,1);
ub = hend*ones(2,1);        

% Initial Guess
 x0 = [1 1]';

% Options
opts = optiset('solver','ipopt','display','iter');

% Build OPTI Problem
Opt = opti('fun',fun,'nl',nlcon,cl,cu,'bounds',lb,ub,'x0',x0,'options',opts);

% Solve NLP
[x,fval,exitflag,info] = solve(Opt);
A_opt(k) = x(1);
C_opt(k) = x(2);
h_opt(k) = h;

Volopt(k) = (1/3).*pi.*h_opt(k).*(A_opt(k)./C_opt(k)).^2.*(3*C_opt(k).^2+h_opt(k).^2);
Areopt(k) = pi*(A_opt(k)./C_opt(k)).^2.*(h.^2+C_opt(k).^2);
rho(k) = m./(pi.*h.*(3.*C_opt(k).^2+h.^2).*((A_opt(k)./C_opt(k)).^2)./3);
Ixopt(k) = rho(k).*pi.*(A_opt(k)/C_opt(k)).^2.*(h.^5./5 + C_opt(k).^2.*h.^3./3) + rho(k).*(pi./4).*(A_opt(k)./C_opt(k)).^4.*(h.^5./5 + 2.*C_opt(k).^2.*h.^3./3 + C_opt(k).^4.*h);

rho = m./(pi.*h.*(3.*C_opt(k).^2+h.^2).*((A_opt(k)./C_opt(k)).^2)./3);
IZ_min(k) = rho.*(pi/2).*(A_opt(k)./C_opt(k)).^4.*(h.^5./5 + (2.*C_opt(k).^2.*h.^3)./3 + C_opt(k).^4.*h)
end

% Omegaz Law without aerodynamics
% omega0 = 1; % [deg/sec]
% Iz10 = IZ_min(1);
% Iz2 = 270;
% Iz1 = IZ_min;
% omega = (Iz10 + Iz2).*omega0./(Iz1 + Iz2);

omega1 = 1; % [deg/sec]
I0 = 270; % For LAnder
Il0= IZ_min(1); % the initial value of Il of the frame
Il = IZ_min; % the present value of I1 for the frame

I1 = I0 + Il0;
I2 = I0 + Il;
omega = (2*I1-I2)*omega1/I1;

figure(1); plot([t0:tstep:tf],IZ_min,'linewidth',4);ylabel('Iz'); xlabel('t'); %axis([0.7 1.19 140 250])
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 16; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
% yticks([140 155 170 185 200 215 230 245 262])
% ax.YLim = [140 262];

figure(2); plot([t0:tstep:tf],A_opt,'linewidth',4);ylabel('a'); xlabel('t'); %axis([hstar 1.19 -Inf Inf])
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 16; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
% yticks([0.6 0.62 0.64 0.66 0.68 0.7 0.72 0.74 0.76 0.78 0.8 0.82])
% ax.YLim = [0.6 0.82];

figure(3); plot([t0:tstep:tf],C_opt,'linewidth',4);ylabel('c'); xlabel('t'); %axis([0.68 1.2 -Inf Inf])
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 16; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
% yticks([0.665 0.7 0.75 0.8 0.85 0.9 0.95 1 1.05 1.1 1.15 1.2])
% ax.YLim = [0.665 1.2];

figure(4); plot([t0:tstep:tf],Volopt,[t0:tstep:tf],1.8202*ones(1,length([t0:tstep:tf])),'linewidth',4); xlabel('t'); ylabel('Volume'); %axis([hstar 1.19 -Inf Inf])
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 16; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
% yticks([1.5 2.00 2.50 3.00 3.50 4.00 4.5 5])
% ax.YLim = [1.5 5];

figure(5); plot([t0:tstep:tf],Areopt,[t0:tstep:tf],4.52389*ones(1,length([t0:tstep:tf])),[t0:tstep:tf],zeros(1,length([t0:tstep:tf])),'linewidth',4); xlabel('t'); ylabel('Area'); %axis([hstar 1.19 -Inf Inf])
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 16; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
% yticks([0 0.5 1 1.5 2 2.5 3 3.5 4 4.52 5])
% ax.YLim = [0 5];

figure(7); plot([t0:tstep:tf],Ixopt,'linewidth',4);ylabel('Ix'); xlabel('t'); axis([t0 tf -Inf Inf])
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 16; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
% ax.YLim = [350 750]; 
ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;

figure(8); hold on % a and c togather
yyaxis left; plot([t0:tstep:tf],A_opt,'k','linewidth',4); %ylabel('a')
% yticks([0.6 0.62 0.64 0.66 0.68 0.7 0.72 0.74 0.76 0.78 0.8 0.82])
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 16; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
% ax.YLim = [0.6 0.82]; 
yyaxis right; plot([t0:tstep:tf],C_opt,':dk','linewidth',4);%ylabel('c')
% yticks([0.665 0.7 0.75 0.8 0.85 0.9 0.95 1 1.05 1.1 1.15 1.2])
% ax.YLim = [0.665 1.2]; 
ax.XAxis.LineWidth = 3;
ax.YAxis(1).LineWidth = 3;
ax.YAxis(2).LineWidth = 3;
grid on; xlabel('t');
legend('a values.....','c values.....')

figure(9); hold on % Iz, Ix togather
yyaxis left; plot([t0:tstep:tf],IZ_min,'k','linewidth',4); %ylabel('a')
% yticks([140 155 170 185 200 215 230 245 262])
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 16; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
% ax.YLim = [140 262];
yyaxis right; plot([t0:tstep:tf],Ixopt,'-.k','linewidth',4);%ylabel('c')
% yticks([350 400 450 500 550 600 650 700 750]); ax.YLim = [0.665 1.2];
ax.XAxis.LineWidth = 3; ax.YAxis(1).LineWidth = 3; ax.YAxis(2).LineWidth = 3;
xlabel('t'); 
% ax.YLim = [350 750];
legend('Iz values.....','Ix values.....')

figure(10); plot([t0:tstep:tf],omega,'k','linewidth',4);ylabel('\omega_z'); xlabel('t'); hold all
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 16; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
% xticks([0.67 0.72 0.77 0.82 0.87 0.92 0.97 1.02 1.08 1.14 1.2])
ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
% yticks([140 155 170 185 200 215 230 245 262])
% ax.YLim = [140 262]; 
% ax.XLim = [0.67 1.2];

% figure(11); plot(IZ_min,omega,'k','linewidth',4);ylabel('\omega_z'); xlabel('IZmin');
% ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 16; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
% % xticks([0.67 0.72 0.77 0.82 0.87 0.92 0.97 1.02 1.08 1.14 1.2])
% ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
% % yticks([140 155 170 185 200 215 230 245 262])
% % ax.YLim = [140 262]; 
% % ax.XLim = [0.67 1.2];

