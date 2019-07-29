function [ F2 ] = Figuras2( f, M, TPSD, r )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here



F2 = figure('name','Datos y Modelo Inicial')
loglog(f,M,'ob','LineWidth',1.5);hold on;
loglog(f,TPSD','color',[0.06,0.7,0.06]','LineWidth',0.5);grid on
title(['Inversi?n del PSD para $r$ =' ' ',num2str(r)],'FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
ylabel('$PSD$','FontSize', 14,'Rotation',90,'interpreter','latex')
legend('M_{Obs}','PSD_{0}'); 




% % % % % % % % % % % GRAFICAS DE LA FUNCION Fun_M.m % % % % % % % % % % % % % % % % 
%Graficado del cociente G0/G1 con y sin suavizado
% figure('name','Funcion M')
% loglog (f,Datos,'b',f,M_smooth,'r','LineWidth',1)
% title('Espectro de Densidad de Potencia observado','FontSize', 12,'interpreter','latex')
% grid on, legend('M','M_{Smooth}')
% xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
% ylabel('$M \left[rk(w) \rigth]$','FontSize', 14,'FontWeight','bold','Rotation',0,'interpreter','latex')

% figure('name','Funcion M')
% loglog(f,Datos,'b'); grid on
% title('Funcion M','FontSize', 12,'interpreter','latex')
% xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
% ylabel('$ M \left[rk(w) \right]$','FontSize', 14,'Rotation',90,'interpreter','latex')

% figure('name','Funcion M suavizada')
% loglog(f,M_smooth,'r')
% title('Funcion M suavizada','FontSize', 12,'interpreter','latex')
% xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
% ylabel('$M \left[rk(w) \rigth]$','FontSize', 14,'FontWeight')
% % % % % % % % % % % FIN GRAFICAS DE LA FUNCION Fun_M.m % % % % % % % % % % % % % 


end

