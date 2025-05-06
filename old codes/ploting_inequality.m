%% plot inequality
% Here we are ploting all values of a,c,Iz,V,Area

clc;
clear;
close all;

Vmin = 1.820219934;
Smax = 4.523893422;

hend = 3*Vmin/Smax;
h = 0.9:0.001:hend;

%% One-Sheeted Hyperboloid:
a = (sqrt(Vmin.*h.*(-Smax.*h + 3.*Vmin).^2./((Smax.*h.^3 + Vmin.*h.^2 + h.*(-Smax.*h + 3.*Vmin).^2).*h.*pi)));
c = (sqrt(h.^2.*(-Smax.*h + 3.*Vmin)./(3.*(Smax.*h + Vmin))));
rho = 576./((1/3).*pi.*h.*(a./c).^2.*(3*c.^2+h.^2));
Iz = rho.*pi/2.*(a./c).^4.*(h.^5./5 + (2.*c.^2.*h.^3)./3 + c.^4.*h);
Ix = rho.*pi.*(a/c).^2.*(h.^5./5 + c.^2.*h.^3./3) + rho.*(pi./4).*(a./c).^4.*(h.^5./5 + 2.*c.^2.*h.^3./3 + c.^4.*h);
Vol = (1/3).*pi.*h.*(a./c).^2.*(3*c.^2+h.^2);
Are = pi*(a./c).^2.*(h.^2-c.^2);

figure(1); plot(h,a,h,zeros(1,length(h)),'linewidth',2); xlabel('h'); ylabel('Parameter A'); grid on; axis([h(1) h(end) -Inf Inf])
figure(2); plot(h,c,h,zeros(1,length(h)),'linewidth',2); xlabel('h'); ylabel('Parameter C'); grid on; axis([h(1) h(end) -Inf Inf])
figure(3); plot(h,Iz,h,zeros(1,length(h)),'linewidth',2); xlabel('h'); ylabel('Moment Iz'); grid on; axis([h(1) h(end) -Inf Inf])
figure(4); plot(h,Ix,h,zeros(1,length(h)),'linewidth',2); xlabel('h'); ylabel('Moment Ix'); grid on; axis([h(1) h(end) -Inf Inf])
figure(5); plot(h,Vol,h,Vmin*ones(1,length(h)),h,zeros(1,length(h)),'linewidth',2); xlabel('h'); ylabel('Volume'); grid on; axis([h(1) h(end) -Inf Inf])
figure(6); plot(h,Are,h,Smax*ones(1,length(h)),h,zeros(1,length(h)),'linewidth',2); xlabel('h'); ylabel('Area'); grid on; axis([h(1) h(end) -Inf Inf])


%% Two-Sheeted Hyperboloid:
h = (3*Vmin/Smax):0.01:2;
alpha = sqrt(9.*Smax.^2.*h.^2 - 30.*Smax.*Vmin.*h + 9.*Vmin.^2);
% a2 = sqrt(6.*Vmin.*Smax.*(Smax.*h-3.*Vmin+alpha).^2)./(pi.*(3.*Smax.*h-3.*Vmin+alpha).*(3.*Smax.*h+3.*Vmin-alpha).^2);
% c2 = (Smax.*h-3.*Vmin+alpha)/(4.*Smax);
a2=sqrt(6.*(Smax.*h - 3.*Vmin + sqrt(9.*Smax.^2*h.^2 - 30.*Smax.*Vmin.*h + 9.*Vmin.^2)).^2.*Vmin.*Smax./(pi*(3.*Smax.*h - 3.*Vmin + sqrt(9.*Smax.^2.*h.^2 - 30.*Smax.*Vmin.*h + 9.*Vmin.^2)).*(3.*Smax.*h + 3.*Vmin - sqrt(9.*Smax.^2.*h.^2 - 30.*Smax.*Vmin.*h + 9.*Vmin.^2)).^2));
c2=(Smax.*h - 3.*Vmin + sqrt(9.*Smax.^2*h.^2 - 30.*Smax.*Vmin.*h + 9.*Vmin.^2))./(4.*Smax);
rho2 = 576./((pi./3).*(c2-h).^2.*(2.*c2+h).*(a2./c2).^2);
IZ2 = rho2.*(pi./2).*(a2./c2).^4.*((h.^5-c2.^5)./5-2.*c2.^2.*(h.^3-c2.^3)./3+c2.^4.*(h-c2));
Ix2 = rho2.*pi.*(a2./c2).^2.*((-c2.^5 + h.^5)./5 - c2.^2.*(-c2.^3 + h.^3)./3) + rho2.*(pi./4).*(a2./c2).^4.*((-c2.^5 + h.^5)./5 - (2.*c2.^2.*(-c2.^3 + h.^3))./3 + c2.^4.*(h - c2));
Vol2 = (pi./3).*(c2-h).^2.*(2.*c2+h).*(a2./c2).^2;
Are2 = pi*(a2./c2).^2.*(h.^2-c2.^2);

figure(6); plot(h,a2,'linewidth',2); xlabel('h'); ylabel('Parameter A'); grid on;
figure(7); plot(h,c2,'linewidth',2); xlabel('h'); ylabel('Parameter C'); grid on;
figure(8); plot(h,IZ2,'linewidth',2); xlabel('h'); ylabel('Moment Iz'); grid on
figure(9); plot(h,Ix2,'linewidth',2); xlabel('h'); ylabel('Moment Ix'); grid on
figure(10); plot(h,Vol2,'linewidth',2); xlabel('h'); ylabel('Volume'); grid on;
figure(11); plot(h,Are2,'linewidth',2); xlabel('h'); ylabel('Area'); grid on


