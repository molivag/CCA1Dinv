clear; close all;  clc

f = (1E-1:0.1:1.25E1);      %Ancho de banda
r = 3;                      %Radio del arreglo circular

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

        nV = (LonReg*NoReg)/W;      %numero de ventanas totales
n_per_wind = W/Dt;                  %datos por ventana
        fs = 1/Dt;                  %
%        f = (0:n_per_wind-1)*(fs/n_per_wind); % Frequency range


Datos = Fun_M( file, NoEst, NoReg, LonReg, Dt, W, Tras, nV, n_per_wind );


% % % % % % % % % % % % % % MODELO DIRECTO % % % % % % % % % % % % % % %  % 

A=115; B=0.9;         %Expresion que define la forma de la curva
Vp_ini = A.*f.^(-B);  %de la velocidad de fase Vp
  
% V0 = 1000;   %m/s  OPRTIMO 1
% Dv = 10;   %m/s
% sigma = .5; %     OPTIMO 0.5
% Vp_ini = V0 + Dv*exp((-f.^2)./sigma)';

OBS= length(f);
PAR= length(Vp_ini);

per=0.0025;
h = Vp_ini*per;

FM = DirectoCCA( f, r, Vp_ini, OBS )'; %transpuesto solo para visualizacion

% % % % % % % % % % % % % % PROCESO DE INVERSION % % % % % % %  % % % % % %
% 
x=[1:length(Vp_ini)];
       %n x %p
Z=zeros(OBS,PAR);

close all
    k=1;
    l=1;
     Vp1 = Vp_ini;
     Vp2 = Vp_ini;
    for j= 1:PAR
     j  ;
         if j==1
            Vp1(l) = Vp_ini(l)+h(l) ; 
            Vp2(l) = Vp_ini(l)-h(l);
            Fper_1=DirectoCCA(f, r, Vp1, OBS)';                      %respuesta f(x+h)
            Fper_1_2h=DirectoCCA(f, r, Vp2, OBS)';                  %respuesta f(x+2h)

           for i=1:OBS
             Z(i,j)=(-3*FM(i) + 4*Fper_1(i) - Fper_1_2h(i))/(2*h(j));
           end
                  
           
           
           
           
        elseif j==(PAR)
            for i=1:OBS
            Z(i,j)=20;
            end
            
            
         else
            j;
            k=k+1;
            Vp_pluss_Perturbada = Vp_ini;
            Vp_minus_Perturbada = Vp_ini;
        
            Vp_pluss_Perturbada(k) = Vp_ini(k)+h(k);
            Vp_minus_Perturbada(k) = Vp_ini(k)-h(k);
            Fper_2= DirectoCCA( f, r, Vp_pluss_Perturbada, OBS )' ;    %FM con f(x+h)
            Bper_2= DirectoCCA( f, r, Vp_minus_Perturbada, OBS )';
                for i=1:OBS    
                    i;
                    Z(i,j)=(Fper_2(i) - Bper_2(i))/ (2*h(j)) ;
                end
         end
    end

    
    
Matriz= Z
Inversa = inv(Z'*Z)
Determinante= det(Z'*Z)


% % % % % % % % % % % GRAFICAS DE LA FUNCION Fun_M.m % % % % % % % % % % % % % % % % 
%Graficado del cociente G0/G1 con y sin suavizado
% figure('name','Funcion M')
% loglog (f,Datos,'b',f,M_smooth,'r','LineWidth',1)
% title('Espectro de Densidad de Potencia observado','FontSize', 12,'interpreter','latex')
% grid on, legend('M','M_{Smooth}')
% xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
% ylabel('$M \left[rk(w) \rigth]$','FontSize', 14,'FontWeight','bold','Rotation',0,'interpreter','latex')

figure('name','Funcion M')
loglog(f,Datos,'b'); grid on
title('Funcion M','FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
ylabel('$ M \left[rk(w) \right]$','FontSize', 14,'Rotation',90,'interpreter','latex')

% figure('name','Funcion M suavizada')
% loglog(f,M_smooth,'r')
% title('Funcion M suavizada','FontSize', 12,'interpreter','latex')
% xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
% ylabel('$M \left[rk(w) \rigth]$','FontSize', 14,'FontWeight')
% % % % % % % % % % % FIN GRAFICAS DE LA FUNCION Fun_M.m % % % % % % % % % % % % % 


% % % % % % % % % % % GRAFICAS DE LA FUNCION DirectoCCA.m % % % % % % % % % % % % %
figure('name','Modelo Directo PSD')
subplot(2,4,[2 3])
semilogx(f,Vp_ini,'k','LineWidth',1);grid on;
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
subplot(2,4,[5 8]);
loglog(f,FM','color',[0.06,0.7,0.06]','LineWidth',0.5)
grid on;	%legend('J_0/J_1')
title(['Modelado directo del cociente espectral de densidad de potencia para $r$ =' ' ',num2str(r)],'FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
ylabel('$\displaystyle \left[\frac{J_{_0}}{J_{_1} }\right]^2 $','FontSize', 14,'Rotation',0,'interpreter','latex')
% % % % % % % % % % % FIN GRAFICAS DE LA FUNCION DirectoCCA.m % % % % % % % % % % %
