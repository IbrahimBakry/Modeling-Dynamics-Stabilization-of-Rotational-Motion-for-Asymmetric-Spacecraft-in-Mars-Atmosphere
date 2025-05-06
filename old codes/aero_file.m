function [Cxfitted, Cyfitted, Cmfitted]=aero_file(A)
%% Aerodynamic File

alpha = [0 5 10 15 20];
Cy = [0 0.15 0.28 0.42 0.53] ;
Cx = [0.24 0.22 0.185 0.13 0.04];
Cm = [0 -0.02 -0.04 -0.05 -0.06];

Cypoly = polyfit(alpha,Cy,2);
Cyfitted = polyval(Cypoly,A);
% figure(1); plot(alpha,Cy,'o'); hold on
% plot([0:0.1:90],Cyfitted); hold off

Cxpoly = polyfit(alpha,Cx,2);
Cxfitted = polyval(Cxpoly,A);
% figure(2); plot(alpha,Cx,'o'); hold on
% plot([0:0.1:90],Cxfitted); hold off

Cmpoly = polyfit(alpha,Cm,2);
Cmfitted = polyval(Cmpoly,A);
% figure(3); plot(alpha,Cm,'o'); hold on
% plot([0:0.1:90],Cmfitted); hold off

end