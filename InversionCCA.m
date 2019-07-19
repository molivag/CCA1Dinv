clear; close all;  clc

% % % % % % % % % % PARAMETROS PARA DATOS OBSERVADOS % % % % % % %  % % % %
file=load ('registros.dat'); %registros.dat antes lamado todas.dat
%NoEst=input('Numero de estaciones en el arreglo circular: ');
NoEst= 20
%NoReg=input('Numero de registros de ruido: ');
NoReg= 18
%LonReg=input('Longitud de registro (seg): ');
LonReg= 65
%Dt=input('Muestreo: ');
Dt= 0.008
%W=input('Tamanio de la ventana (seg): ');
W= 1
%Tras=input('Traslape de ventanas 1(0%) o 2(50%): ');
Tras= 1

      r = 15;                    %Radio del arreglo circular
     nV = (LonReg*NoReg)/W;      %numero de ventanas totales
n_x_ven = W/Dt                   %datos por ventana             RENOMBRAR POR Nven al parecer siempre ser? igual a las observaciones por lo que se puede quitar y ahorrar una variable
     fs = 1/Dt;                  %Frecuencia Maxima
     ds = 1/(n_x_ven*Dt);        %Muestreo frec
      f = (0:ds:n_x_ven-1)*(fs/n_x_ven); % Frequency range
%      f = (1E-1:0.9:1.25E1);    %Ancho de banda

Datos = Fun_M( file, NoEst, NoReg, LonReg, Dt, W, Tras, nV, n_x_ven );


% % % % % % % % % % % % % % MODELO DIRECTO % % % % % % % % % % % % % % %  % 

A=1000; B=0.9;          %Expresion que define la forma de la curva
Vp = A.*f.^(-B);        %de la velocidad de fase Vp
disp('Tamanio de la matriz')
OBS = length(Datos)     %n
PAR = length(Vp);       %p
Z=zeros(OBS,PAR);
dim=size(Z)
Elementos=OBS*PAR
per=0.0025;             %Perturbacion en el Jacobiano
FM = DirectoCCA( f, r, Vp, OBS )'; %transpuesto solo para visualizacion

% % % % % % % % % % % % % % MODELO INVERSION % % % % % % %  % % % % % %









disp('Verificaci?n de la estabilidad de la matriz ZtZ')
INV_ZtZ=inv(Z'*Z);
Determinante=det(Z'*Z)
Norma=norm(INV_ZtZ)




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





% % % % % % % % % % % GRAFICAS DEL JACOBIANO % % % % % % % % % % % % % % % % 

% figure('name','Modelo Directo')
% loglog(f,Resp,'k'); grid on; title('Modelo Directo')
%
% figure('name','Jacobisno')
% plot(Z); title('Matriz de Sensibilidad Z')




% % % % % % % % % % % GRAFICAS DE LA FUNCION DirectoCCA.m % % % % % % % % % % % % %
figure('name','Modelo Directo PSD')
subplot(3,4,[2 3])
semilogx(f,Vp,'k','LineWidth',1);grid on;
title('Velocidad de fase $\nu_p$ inicial','FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
ylabel('$\nu_p$','FontSize', 14,'FontWeight','bold','Rotation',0,'interpreter','latex')
% subplot(3,1,2);
% loglog(f,J0kr,'r',f,J1kr,'b','LineWidth',0.5)
% title('Funciones de Bessel de primera especie para  $\nu \in [0, 1]$','FontSize', 12,'interpreter','latex')
% xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
% ylabel('$J_\nu[rk(\omega)]$','FontSize', 15,'FontWeight','bold','Rotation',90,'interpreter','latex')
% grid on; etiqueta = {'$$ \nu = 0 $$','$$ \nu = 1$$'};
% legend(etiqueta, 'Interpreter','latex', 'Location','best')%,'FontSize', 14)
subplot(3,4,[5 8]);
loglog(f,Datos(1:length(f)),'b'); grid on
title('Funcion M: Datos Observados','FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
ylabel('$ M \left[rk(w) \right]$','FontSize', 14,'Rotation',90,'interpreter','latex')

subplot(3,4,[9 12]);
loglog(f,FM','color',[0.06,0.7,0.06]','LineWidth',0.5)
%loglog(f,FM','color',[0.6,0.1,0]','LineWidth',0.5)
grid on;	%legend('J_0/J_1')
title(['Modelado directo del cociente espectral de densidad de potencia para $r$ =' ' ',num2str(r)],'FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
ylabel('$\displaystyle \left[\frac{J_{_0}}{J_{_1} }\right]^2 $','FontSize', 14,'Rotation',0,'interpreter','latex')
% % % % % % % % % % % FIN GRAFICAS DE LA FUNCION DirectoCCA.m % % % % % % % % % % %
