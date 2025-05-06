%% this program find the minimum values of Iz1 for range h a c

clc;
clear;
close all;

Vmin = 1.820219934;
Smax = 4.523893422;
hend = 3*Vmin/Smax;
hstart = 0.7;
hstep = 0.01;

k = 0;
for h = hstart:hstep:hend
    k = k + 1;
    a = (sqrt(Vmin.*h.*(-Smax.*h + 3.*Vmin).^2./((Smax.*h.^3 + Vmin.*h.^2 + h.*(-Smax.*h + 3.*Vmin).^2).*h.*pi)));
    c = (sqrt(h.^2.*(-Smax.*h + 3.*Vmin)./(3.*(Smax.*h + Vmin))))
    ii = 0;
    for i = a:(a/50):5*a
        ii = ii + 1;
        jj = 0;
        for j = c:(c/50):5*c
            jj = jj + 1;
            IZ(ii,jj) = (576./(pi./3*(2.*c + h).*(c - h).^2.*(a./c).^2)).*pi/2.*(a./c).^4.*(h.^5./5 + (2.*c.^2.*h.^3)./3 + c.^4.*h);
        end
    end
    IZZ(k) = min(IZ(:));
    [Row, Column] = ind2sub( size(IZ), find(IZ == IZZ(k), 1 ));
    AA = a:(a/50):5*a;
    CC = c:(c/50):5*c;
    A(k) = AA(Row);
    C(k) = CC(Column);
end

% Checking Volume and Area
h = hstart:hstep:hend;
Vol = (1/3).*pi.*h.*(A./C).^2.*(3*C.^2+h.^2);
Are = pi*(A./C).^2.*(h.^2-C.^2);

figure(1); plot(h,IZZ,'linewidth',2);grid on;ylabel('IZ'); xlabel('h')
figure(2); plot(h,A,'linewidth',2);grid on;ylabel('A'); xlabel('h')
figure(3); plot(h,C,'linewidth',2);grid on;ylabel('C'); xlabel('h')
figure(4); plot(h,Vol,h,1.8202*ones(1,length(h)),'linewidth',2); grid on; xlabel('h'); ylabel('Volume')
figure(5); plot(h,Are,h,4.52389*ones(1,length(h)),'linewidth',2); grid on; xlabel('h'); ylabel('Area')
figure(6); plot(A,C,'linewidth',2); grid on; xlabel('A'); ylabel('C')

