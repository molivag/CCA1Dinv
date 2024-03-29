function [ F7 ] = Fig7( finv, Xmc )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

F7 =figure(4); 
set( F7,'Position', [700 10 450 310],'name','Velocidad de fase calculada');
plot(finv,Xmc,'+r','LineWidth',1.5)
title('Velocidad de fase calculada')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
y=ylabel('$\nu_p$','FontSize', 12,'FontWeight','bold','Rotation',0,'interpreter','latex');
set(y, 'Units', 'Normalized', 'Position', [-0.085, 0.5, 0])
grid on
end

