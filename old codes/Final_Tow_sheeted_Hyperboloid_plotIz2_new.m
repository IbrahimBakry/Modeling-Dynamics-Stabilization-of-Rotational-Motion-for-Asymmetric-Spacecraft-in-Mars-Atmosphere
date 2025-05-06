%% this program find the minimum values of Iz1 for range h a c
clc; clear; close all;

Vmin = 1.820219934;
Smax = 4.523893422;
hstar= 1.98;%1.26
hend = 2.5;
hstep= 0.06;

% Parameters Intervals
a = [0.01:0.005:5];
c = [0.01:0.005:5];

k = 0; 
for h = hstar:hstep:hend
    k = k + 1; ii = 0; q = 0;
    for i = 1:length(a)
        ii = ii + 1; jj = 0;
        for j = 1:length(c)
            jj = jj + 1;          
            Vol = (pi./3).*(c(j)-h).^2.*(2.*c(j)+h).*(a(i)./c(j)).^2;
            Are = pi*(a(i)./c(j)).^2.*(h.^2+c(j).^2);

            % Checking the condition and selecting a & c
            if (Vol >= Vmin) && (Are <= Smax) && (Are > 0)
                q = q + 1;
                A(q) = a(i);
                C(q) = c(j);
                IZ(q) = (576./((pi./3).*(C(q)-h).^2.*(2.*C(q)+h).*(A(q)./C(q)).^2)).*(pi./2).*(A(q)./C(q)).^4.*((h.^5-C(q).^5)./5-2.*C(q).^2.*(h.^3-C(q).^3)./3+C(q).^4.*(h-C(q)));
                IZAC(q,:) = [A(q), C(q), IZ(q)];
            end
        end
    end
    if q <= 1
        disp('There are no values of a & c that satisfy the Conditions')
    else
    % Getting the minimum value
    IZ_min(k) = min(IZAC(:,3));
    IZ_max(k) = max(IZAC(:,3));
   
    val = find(IZ==IZ_min(k));
    A_opt(k) = A(val);
    C_opt(k) = C(val);
    h(k)=h;
    end
end

Vol = (pi./3).*(C_opt-h).^2.*(2.*C_opt+h).*(A_opt./C_opt).^2;
Are = pi*(A_opt./C_opt).^2.*(h.^2+C_opt.^2);
c2=C_opt; a2=A_opt;
rho2 = 576./((pi./3).*(c2-h).^2.*(2.*c2+h).*(a2./c2).^2);
Ix2 = rho2.*pi.*(a2./c2).^2.*((-c2.^5 + h.^5)./5 - c2.^2.*(-c2.^3 + h.^3)./3) + rho2.*(pi./4).*(a2./c2).^4.*((-c2.^5 + h.^5)./5 - (2.*c2.^2.*(-c2.^3 + h.^3))./3 + c2.^4.*(h - c2));

figure(1); plot([hstar:hstep:hend],IZ_min,'linewidth',4);ylabel('Iz'); xlabel('h'); %axis([0.7 1.19 140 250])
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 16; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
xticks([1.98 2.04 2.1 2.16 2.22 2.28 2.34 2.4 2.46])
ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
% yticks([140 155 170 185 200 215 230 245 262])
% ax.YLim = [140 262]; 
ax.XLim = [1.98 2.46];

figure(2); plot([hstar:hstep:hend],A_opt,'linewidth',4);ylabel('a'); xlabel('h'); %axis([hstar 1.19 -Inf Inf])
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 16; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
% xticks([0.67 0.72 0.77 0.82 0.87 0.92 0.97 1.02 1.08 1.14 1.2])
ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
% yticks([0.6 0.62 0.64 0.66 0.68 0.7 0.72 0.74 0.76 0.78 0.8 0.82])
% ax.YLim = [0.6 0.82]; 
ax.XLim = [1.98 2.46];

figure(3); plot([hstar:hstep:hend],C_opt,'linewidth',4);ylabel('c'); xlabel('h'); %axis([0.68 1.2 -Inf Inf])
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 16; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
xticks([1.98 2.04 2.1 2.16 2.22 2.28 2.34 2.4 2.46])
ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
% yticks([0.665 0.7 0.75 0.8 0.85 0.9 0.95 1 1.05 1.1 1.15 1.2])
% ax.YLim = [0.665 1.2]; 
ax.XLim = [1.98 2.46];

figure(4); plot([hstar:hstep:hend],Vol,[hstar:0.01:hend],1.8202*ones(1,length([hstar:0.01:hend])),'linewidth',4); xlabel('h'); ylabel('Volume'); %axis([hstar 1.19 -Inf Inf])
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 16; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
xticks([1.98 2.04 2.1 2.16 2.22 2.28 2.34 2.4 2.46])
ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
ax.YLim = [1.5 3.5]; 
ax.XLim = [1.98 2.46];
yticks([1.5 1.75 1.82 2.00 2.25 2.50 2.75 3.00 3.25 3.5])

figure(5); plot([hstar:hstep:hend],Are,[hstar:0.01:hend],4.52389*ones(1,length([hstar:0.01:hend])),[hstar:0.01:hend],zeros(1,length([hstar:0.01:hend])),'linewidth',4); xlabel('h'); ylabel('Area'); %axis([hstar 1.19 -Inf Inf])
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 16; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
xticks([1.98 2.04 2.1 2.16 2.22 2.28 2.34 2.4 2.46])
ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
% yticks([0 0.5 1 1.5 2 2.5 3 3.5 4 4.52])
ax.YLim = [0 4.524]; 
ax.XLim = [1.98 2.46001];

figure(7); plot([hstar:hstep:hend],Ix2,'linewidth',4);ylabel('Ix'); xlabel('h'); %axis([hstar 1.19 -Inf Inf])
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 16; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
% ax.YLim = [350 750]; ax.XLim = [0.67 1.2];
xticks([1.98 2.04 2.1 2.16 2.22 2.28 2.34 2.4 2.46])
ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
ax.XLim = [1.98 2.46];

figure(8); % a and c togather
yyaxis left; plot([hstar:hstep:hend],A_opt,'linewidth',4); %ylabel('a')
% yticks([0.6 0.62 0.64 0.66 0.68 0.7 0.72 0.74 0.76 0.78 0.8 0.82])
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 16; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
% ax.YLim = [0.6 0.82]; 
yyaxis right; plot([hstar:hstep:hend],C_opt,':d','linewidth',4);%ylabel('c')
% yticks([0.665 0.7 0.75 0.8 0.85 0.9 0.95 1 1.05 1.1 1.15 1.2])
% ax.YLim = [0.665 1.2]; 
ax.XAxis.LineWidth = 3;
ax.YAxis(1).LineWidth = 3;
ax.YAxis(2).LineWidth = 3;
xlabel('h');
xticks([1.98 2.04 2.1 2.16 2.22 2.28 2.34 2.4 2.46])
legend('a values.....','c values.....')
ax.XLim = [1.98 2.46];

figure(9); % Iz, Ix togather
yyaxis left; plot([hstar:hstep:hend],IZ_min,'linewidth',4); %ylabel('a')
% yticks([140 155 170 185 200 215 230 245 262])
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.7; ax.FontSize = 16; ax.FontWeight= 'bold'; ax.LineWidth = 0.8; 
% ax.YLim = [140 262];
yyaxis right; plot([hstar:hstep:hend],Ix2,'-.','linewidth',4);%ylabel('c')
% yticks([350 400 450 500 550 600 650 700 750]); ax.YLim = [0.665 1.2];
ax.XAxis.LineWidth = 3; ax.YAxis(1).LineWidth = 3; ax.YAxis(2).LineWidth = 3;
xlabel('h'); 
% ax.YLim = [350 750]; 
ax.XLim = [1.98 2.46];
xticks([1.98 2.04 2.1 2.16 2.22 2.28 2.34 2.4 2.46])
legend('Iz values.....','Ix values.....')

