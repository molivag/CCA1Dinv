function [ F1] = Figuras( f, Vp, M, TPSD, r )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here



% % % % % % % % % % % GRAFICAS DE LA FUNCION DirectoCCA.m % % % % % % % % % % % % %
F1= figure('name','Modelo Directo PSD');
subplot(3,4,[2 3])
semilogx(f,Vp,'k','LineWidth',1);grid on;
title('Velocidad de fase $\nu_p$ inicial','FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
ylabel('$\nu_p$','FontSize', 14,'FontWeight','bold','Rotation',0,'interpreter','latex')

subplot(3,4,[5 8]);
loglog(f,M(1:length(f)),'b'); grid on
%loglog(f,M','b'); grid on
title('Funcion M: Datos Observados','FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
ylabel('$ M \left[rk(w) \right]$','FontSize', 14,'Rotation',90,'interpreter','latex')

subplot(3,4,[9 12]);
loglog(f,TPSD(1:length(f))','color',[0.06,0.7,0.06]','LineWidth',0.5)
%loglog(f,TPSD','color',[0.06,0.7,0.06]','LineWidth',0.5)
title(['Modelado directo del cociente espectral de densidad de potencia para $r$ =' ' ',num2str(r)],'FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex');grid on
ylabel('$\displaystyle \left[\frac{J_{_0}}{J_{_1} }\right]^2 $','FontSize', 14,'Rotation',0,'interpreter','latex')

end

