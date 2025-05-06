%% Charting hyperboloid
close all;

Rbase = Rmax;
ybase = -Rbase:0.01:Rbase;
l1 = zeros(size(ybase));

figure(20); 
plot(ybase,l1,'k','linewidth',2); hold on % Charting the first base line 

for i = 1 : length(A_opt)
    a = A_opt(i);
    c = C_opt(i);
    y = a:0.0001:Rbase;
    
    z10 = (c./a).*sqrt(y.^2 - a^2)/2;
    z1 = z10(end) - (c./a).*sqrt(y.^2 - a^2)/2;
    l1 = z1(end)*ones(size(-Rbase:0.01:Rbase));
    l2 = z10(end)*ones(size(-a:0.001:a));
    
%     if z10(end) >= hmin && z10(end) <= hmax
       plot(y,z1,'-k',-y,z1,'-k','linewidth',0.5); grid on; hold all; % charting the side curves (parabola)
       plot(-Rbase:0.01:Rbase,l1,'-k','linewidth',0.5); hold all; % charting the side curves (parabola)
       plot(-a:0.001:a,l2,'-k','linewidth',0.5); hold all; % charting the side curves (parabola)
%         axis([-1.5 1.5 -2 0])
       pause(0.5)
%     end
end

