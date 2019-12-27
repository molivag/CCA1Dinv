function [ F1] = Fig1(f, M )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here



% % % % % % % % % % % GRAFICAS DE LA FUNCION DirectoCCA.m % % % % % % % % % % % % %
figure('name','Datos observados y Modelo Inicial')
F1 = subplot(3,7,1:7);
 loglog(f,M,'b'); grid on
title('Funcion M: Datos Observados','FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
y=ylabel('$ M \left[rk(w) \right]$','FontSize', 11,'Rotation',90,'interpreter','latex');
set(y, 'Units', 'Normalized', 'Position', [-0.04, 0.5, 0])

end
