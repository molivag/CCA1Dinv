clear; %close all

f = (1E-1:0.1:10E1);            %Rango de frecuencia
r = 3; %*ones(1,length(w));    %Radio del arreglo


%A=4; B=0.08;      %Expresion que define la forma de la curva
%Vp = A.*f.^(-B);  %de la velocidad de fase Vp
  
V0 = 100;   %m/s  OPRTIMO 1
Dv = 10;   %m/s
sigma = .5; %     OPTIMO 0.5
Vp = V0 + Dv*exp((-f.^2)./sigma);



kr = r*(2*pi*f./Vp);      %Numero de onda K   

J0kr = (besselj(0,kr)).^2;
J1kr = (besselj(1,kr)).^2;

FMkr = J0kr./J1kr;  %Forward Model


figure('name','Modelo Directo PSD')
%figure('name','Velocidad de fase inicial')
subplot(3,1,1);
semilogx(f,Vp,'k','LineWidth',1);grid on;
title('Velocidad de fase $\nu_p$ inicial','FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
ylabel('$\nu_p$','FontSize', 14,'FontWeight','bold','Rotation',0,'interpreter','latex')

subplot(3,1,2);
loglog(f,J0kr,'r',f,J1kr,'b','LineWidth',0.5)
title('Funciones de Bessel de primera especie para  $\nu \in [0, 1]$','FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
ylabel('$J_\nu[rk(\omega)]$','FontSize', 15,'FontWeight','bold','Rotation',90,'interpreter','latex')
grid on; etiqueta = {'$$ \nu = 0 $$','$$ \nu = 1$$'};
legend(etiqueta, 'Interpreter','latex', 'Location','best')%,'FontSize', 14)

subplot(3,1,3);
loglog(f,FMkr,'color',[0.06,0.7,0.06]','LineWidth',0.5)
grid on;	%legend('J_0/J_1')
title(['Modelado directo del cociente espectral de densidad de potencia para $r$ =' ' ',num2str(r)],'FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
ylabel('$\displaystyle \left[\frac{J_{_0}}{J_{_1} }\right]^2 $','FontSize', 14,'Rotation',0,'interpreter','latex')