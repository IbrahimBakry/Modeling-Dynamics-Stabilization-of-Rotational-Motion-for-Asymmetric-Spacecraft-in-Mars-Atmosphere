%% Comparing Volumes

hopt = 0.7;
A_opt = 0.5;
C_opt = 0.5;
Volopt1 = (1/3).*pi.*hopt.*(A_opt./C_opt).^2.*(3*C_opt.^2+hopt.^2)
Volopt2 = (pi.*A_opt.*C_opt.*log((C_opt.^2).^(1/2)))./2 + (A_opt.*pi.*((hopt.*(C_opt.^2 + hopt.^2).^(1/2))./2 - (C_opt.^2.*log((C_opt.^2 + hopt.^2).^(1/2) - hopt))./2))/C_opt

