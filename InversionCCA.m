clear; close all;  clc

f = (1E-1:0.1:1E1);            %Rango de frecuencia
r = 3; %*ones(1,length(w));   %Radio del arreglo

% 
% file=load ('registros.dat');
% NoEst=input('Numero de estaciones en el arreglo circular: ');
% NoReg=input('Numero de registros de ruido: ');
% LonReg=input('Longitud de registro (seg): ');
% Dt=input('Muestreo: ');
% W=input('Tamanio de la ventana (seg): ');
% Trsl=input('Traslape de ventanas 1(0%) o 2(50%): ');
% 


A=115; B=0.9;      %Expresion que define la forma de la curva
Vp_ini = A.*f.^(-B);  %de la velocidad de fase Vp
  
% V0 = 1000;   %m/s  OPRTIMO 1
% Dv = 10;   %m/s
% sigma = .5; %     OPTIMO 0.5
% Vp_ini = V0 + Dv*exp((-f.^2)./sigma)';

OBS= length(f);
PAR= length(Vp_ini);

per=0.0025;
h = Vp_ini*per;

% % % % % % % % % % M O D E L O   D I R E C T O % % % % % % %  % % % % % %
FM = DirectoCCA( f, r, Vp_ini, OBS )'; %transpuesto solo para visualizacion


% % % % % % % % % % I N V E R S I O N  I T E R A T I V A % % % % % % %  % %

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

    
    




% % % % % % % % % % % % % % GRAFICAS DEL PROGRAMA FUN M % % % % % % % % % % % % % % % % 
%Graficado del cociente G0/G1 con y sin suavizado
figure('name','Funcion M')
loglog (f,M,'b',f,M_smooth,'r');title('M[rk(w)]')
legend('M','M_{Smooth}')
grid on; xlabel('f'); ylabel('M')

figure('name','Funcion M')
loglog(f,M,'b'); grid on

figure('name','Funcion M suavizada')
loglog(f,M_smooth,'r')

% % % % % % % % % % % % % % FIN GRAFICAS DEL PROGRAMA FUN M % % % % % % % % % % % % % %





% % % % % % % % % % % % % F I G U R A S % % % % % % %  % % % % % % % % % % 
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

