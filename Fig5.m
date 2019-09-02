function [ F5 ] = Fig5( finv, M2, TPSDR, r )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here



F5 =figure(2);
set( F5, 'Position', [150 400 1000 300], 'Name', 'Ajustes');
loglog(finv,M2,'ob','LineWidth',1.5);hold on;
loglog(finv,TPSDR','color',[0.06,0.7,0.06]','LineWidth',0.5);grid on
title(['Inversion del PSD para un radio =' ' ',num2str(r),' ' 'm'],'FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
ylabel('$PSD$','FontSize', 14,'Rotation',90,'interpreter','latex')
legend('M_{Obs}','PSD_{0}', 'Location', 'Best'); 



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

