
poly_A = polyfit([hstar:hstep:hend],A_opt,2);
ya = polyval(poly_A,[hstar:hstep:hend]);
figure(1); plot([hstar:hstep:hend],A_opt,'o'); hold on
plot([hstar:hstep:hend],ya); hold off

poly_C = polyfit([hstar:hstep:hend],C_opt,2);
yc = polyval(poly_C,[hstar:hstep:hend]);
figure(2); plot([hstar:hstep:hend],C_opt,'o'); hold on
plot([hstar:hstep:hend],yc); hold off

poly_Iz = polyfit([hstar:hstep:hend],IZ_min,2);
yiz = polyval(poly_Iz,[hstar:hstep:hend]);
figure(3); plot([hstar:hstep:hend],IZ_min,'o'); hold on
plot([hstar:hstep:hend],yiz); hold off

poly_Ix = polyfit([hstar:hstep:hend],Ix,3);
yix = polyval(poly_Ix,[hstar:hstep:hend]);
figure(4); plot([hstar:hstep:hend],Ix,'o'); hold on
plot([hstar:hstep:hend],yix); hold off

