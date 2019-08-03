function [ F8 ] = Fig8( finv, Vp, TPSD, r, F1, F3)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


hold on
delete(F3)
F8 = subplot(4,4,[10 11]);
semilogx(finv,Vp,'LineWidth',1);grid on;
title('Velocidad de fase $\nu_p$ inicial','FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
y=ylabel('$\nu_p$','FontSize', 13,'FontWeight','bold','Rotation',0,'interpreter','latex');
set(y, 'Units', 'Normalized', 'Position', [-0.18, 0.4, 0])
%linkaxes([F3, F8],'x')


F9 = subplot(4,4,[13 16]);
hold on
delete(F9)
F9 = subplot(4,4,[13 16]);
loglog(finv,TPSD','color',[0.06,0.7,0.06]','LineWidth',0.5)
title(['Modelado directo del cociente espectral de densidad de potencia para $r$ =' ' ',num2str(r)],'FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex');grid on
y=ylabel('$\displaystyle \left[\frac{J_{_0}}{J_{_1} }\right]^2 $','FontSize', 11,'Rotation',0,'interpreter','latex');
set(y, 'Units', 'Normalized', 'Position', [-0.08, 0.3, 0])
linkaxes([F1,F9],'x')

end

