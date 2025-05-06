%% this program find the minimum values of Iz2 for range h a c

clc;
clear;
close all;

Vmin = 1.820219934;
Smax = 4.523893422;
hend = 3*Vmin/Smax;
k = 0;
for h = hend:0.01:50
    k = k + 1;
    a = sqrt(6.*(Smax.*h - 3.*Vmin + sqrt(9.*Smax.^2*h.^2 - 30.*Smax.*Vmin.*h + 9.*Vmin.^2)).^2.*Vmin.*Smax./(pi*(3.*Smax.*h - 3.*Vmin + sqrt(9.*Smax.^2.*h.^2 - 30.*Smax.*Vmin.*h + 9.*Vmin.^2)).*(3.*Smax.*h + 3.*Vmin - sqrt(9.*Smax.^2.*h.^2 - 30.*Smax.*Vmin.*h + 9.*Vmin.^2)).^2));
    c = (Smax.*h - 3.*Vmin + sqrt(9.*Smax.^2*h.^2 - 30.*Smax.*Vmin.*h + 9.*Vmin.^2))./(4.*Smax);
    ii = 0;
    for i = 0:(a/50):a
        ii = ii + 1;
        jj = 0;
        for j = 0:(c/50):c
            jj = jj + 1;
            IZ(ii,jj) = (576./((pi./3).*(c-h).^2.*(2.*c+h).*(a./c).^2)).*(pi./2).*(a./c).^4.*((h.^5-c.^5)./5-2.*c.^2.*(h.^3-c.^3)./3+c.^4.*(h-c));
        end
    end
    IZZ(k) = max(IZ(:));
    [Row, Column] = ind2sub( size(IZ), find(IZ == IZZ(k), 1 ));
    AA = 0:(a/50):a;
    CC = 0:(c/50):c;
    A(k) = AA(Row);
    C(k) = CC(Column);
end

% Checking Volume and Area
h = hend:0.01:50;
Vol = (pi./3).*(c-h).^2.*(2.*c+h).*(a./c).^2;
Are = pi*(a./c).^2.*(h.^2-c.^2);

figure(1); plot(h,IZZ,'linewidth',2);grid on;ylabel('IZ'); xlabel('h')
figure(2); plot(h,A,'linewidth',2);grid on;ylabel('A'); xlabel('h')
figure(3); plot(h,C,'linewidth',2);grid on;ylabel('C'); xlabel('h')
figure(4); plot(h,Vol,h,1.8202*ones(1,length(h)),'linewidth',2); grid on; xlabel('h'); ylabel('Volume')
figure(5); plot(h,Are,h,4.52389*ones(1,length(h)),'linewidth',2); grid on; xlabel('h'); ylabel('Area')
figure(6); plot(A,C,'linewidth',2); grid on; xlabel('A'); ylabel('C')


% %% 
% h = 0.01:0.1:hend;
% a = (sqrt(Vmin.*h.*(-Smax.*h + 3.*Vmin).^2./((Smax.*h.^3 + Vmin.*h.^2 + h.*(-Smax.*h + 3.*Vmin).^2).*h.*pi)));
% c = (sqrt(h.^2.*(-Smax.*h + 3.*Vmin)./(3.*(Smax.*h + Vmin))));
% figure(7); plot(h,a)
% figure(8); plot(h,c)