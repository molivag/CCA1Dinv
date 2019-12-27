function [ F3, Opcion, anss] = Fig3( finv, Vp, TPSDR, r, F1, F2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

F3 = subplot(3,7,[12:14 19:21]);
semilogx(finv,Vp,'+k','LineWidth',1);grid on;
title('Velocidad de fase $\nu_p$ inicial','FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
y=ylabel('$\nu_p$','FontSize', 13,'FontWeight','bold','Rotation',0,'interpreter','latex');
set(gca,'YAxisLocation','right');
set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0])

% linkaxes([F1,F2,F3],'x')


F4 = subplot(3,7,15:18);
loglog(finv,TPSDR','color',[0.06,0.7,0.06]','LineWidth',0.5)
% title(['Modelo directo del cociente del espectro de densidad de potencia con radio =' ' ',num2str(r), ' ' ' m'], 'FontSize', 12,'interpreter','latex')
title(['Modelo directo de la funcion M con radio =' ' ',num2str(r), ' ' ' m'], 'FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex');grid on
y=ylabel('$\displaystyle \left[\frac{J_{_0}}{J_{_1} }\right]^2 $','FontSize', 11,'Rotation',0,'interpreter','latex');
set(y, 'Units', 'Normalized', 'Position', [-0.11, 0.4, 0])
%set(y, 'Units', 'Normalized', 'Position', [-0.08, 0.3, 0])
% set(gca,'YAxisLocation','right');
linkaxes([F1,F2,F4],'xy')



Opcion = input([' + + + + + + + Tipo de Inversion + + + + + + + +'...
'\n +                                             + ',...
'\n + - - - 1.- Minimo Cuadrado Estandar: - - - - + ',...
'\n + - - - 2.- Regularizacion de Tikhonov: - - - + ',...
'\n +                                             + ',...
'\n + + + + + + + + + + + + + + + + + + + + + + + + ',...
'\n Opcion = ']);
if Opcion == 1
    disp ('Inversion por Minimo Cuadrado Estandar')
elseif Opcion == 2
    disp('Inversi?n por Regularizacion de Tikhonov')
end

  anss = questdlg('      Proceder con la Inversion?', ...
'Proceso Completado','Yes','No','No');

end

