function [ F8 ] = Fig8( finv, Vp, TPSDR, r, F1)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

F8 = subplot(3,7,[12:14 19:21]);
hold on
delete(F8)
F8 = subplot(3,7,[12:14 19:21]);
semilogx(finv,Vp,'+k','LineWidth',1);grid on;
title('Velocidad de fase $\nu_p$ inicial','FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
y=ylabel('$\nu_p$','FontSize', 13,'FontWeight','bold','Rotation',0,'interpreter','latex');
set(gca,'YAxisLocation','right');
set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0])
%linkaxes([F3, F8],'x')


F9 = subplot(3,7,15:18);
hold on
delete(F9)
F9 = subplot(3,7,15:18);
loglog(finv,TPSDR','color',[0.06,0.7,0.06]','LineWidth',0.5)
title(['Modelo directo de la funcion M con radio =' ' ',num2str(r), ' ' ' m'], 'FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex');grid on
y=ylabel('$\displaystyle \left[\frac{J_{_0}}{J_{_1} }\right]^2 $','FontSize', 11,'Rotation',0,'interpreter','latex');
set(y, 'Units', 'Normalized', 'Position', [-0.11, 0.4, 0])
linkaxes([F1,F9],'xy')

end

