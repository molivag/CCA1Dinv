function [ F2 ] = Fig2( finv, M2, F1 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


F2 = subplot(3,7,8:11);
loglog(finv,M2,'r'); grid on
title('Ancho de banda de datos observadoa invertir','FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
y=ylabel('$ M \left[rk(w) \right]$','FontSize', 11,'Rotation',90,'interpreter','latex') ;
set(y, 'Units', 'Normalized', 'Position', [-0.08, 0.5, 0])
linkaxes([F1,F2],'xy') 
end

