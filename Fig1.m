function [ F1] = Fig1ls( f, M )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here



% % % % % % % % % % % GRAFICAS DE LA FUNCION DirectoCCA.m % % % % % % % % % % % % %
figure('name','no se')
F1 = subplot(4,4,[1 4]);
loglog(f(1:length(M)),M,'b'); grid on
title('Funcion M: Datos Observados','FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
ylabel('$ M \left[rk(w) \right]$','FontSize', 14,'Rotation',90,'interpreter','latex')

end

