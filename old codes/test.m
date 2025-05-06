%test
clc;clear;

k = 0;
A = 0.66:-0.01:0.54;
for h_opt = 1.5:0.04:2
    k = k + 1; hh(k) = h_opt;
    C_opt = 2;
    A_opt = A(k);
Volopt(k) = (1/3).*pi.*h_opt.*(A_opt./C_opt).^2.*(3*C_opt.^2+h_opt.^2);
end
figure; plot(hh,Volopt)