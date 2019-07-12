  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                     *** TEOR??A DE INVERSI??N DE DATOS ***                       %   
  %                                                                                %
  %      Proyecto de clase: Soluci???n al problema inverso del SEV                 %      
  %                                                                                %
  %  Programa de inversi??n para el problema de la inyecci??n de corriente directa   %
  %  en medios estratificados 1-D.                                                 %
  %                                                                                %
  %  El codigo invierte unicamente valores de resistividad y esta limitado para    % 
  %  un modelo de n<=5 capas, de ser necesarias menos de 5, deben modificarse la   %
  %  funci???n Jacobiano.                                                          %
  %                                                                                %
  %  Linversi???n se realiza empleando 4 esquemas de estimaci???n:                 %
  %                                                                                %   
  %     1) M???nimo cuadrado estandar.                                             %
  %     2) M???nimo cuadrado ponderado.                                            %
  %     3) M???nimo cuadrado regularizado.                                         %
  %     4) Descomposici???n en valores singulares.                                 %
  %                                                                                %
  %  El c???digo realiza  3 itraciones en los estimadores estandar y ponderado, y 5%
  %  iteraciones en el estimador recursivo. Reaaliza adem???s la DVS  de la matriz %
  %  de Sensibilidad J. En la iteracion final de cada estimador, se calcula  la    %
  %  nueva matriz J a partir de la respuesta al modelo estimado Y(Xmc) y al modelo % 
  %  mismo (Xmc). El prblema directo se resuelve en cada iteraci???n.              %
  %                                                                                %
  %  Finalmente se generan 4 figuras en las que se muestran los mejores            %
  %  ajustes, modelos, residuales y Jacobianos de las ultimas iteraciones          %
  %  de cada estimador.                                                            %
  %                                                                                %
  %  Los resultados de la inversi???n estan fuertemente influenciados por el       %
  %  modelo inicial "rho" y se recomienda optimizarse empleando el Log. de         %
  %  la resistividad y esquemas definidos para problemas mal condicionados.        %  
  %                                                                                %
  %                                                                                %
  %                                                                                %
  %  Marco Antonio Oliva Guti??rrez                           Geof??sica Aplicada    % 
  %                                                                                %
  %                                                                                %
  %  V 7.5                                                           08/08/2017    %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear all, close all
%Numero de mediciones
m = 17;

%Medciones en AB/2
ab=[ 3 5 7 10 13 16 20 25 32 40 50 65 80 100 130 160 200];
sev7=[29.47 23.75 19.44 14.70 13.07 11.00 9.10 8.02...
      7.25 7.05 8.83 9.25 11.06 13.57 17.54 20.47 27.19]'; %sinteticos

for ii=1:m
sev7_ruido(ii)= sev7(ii)*(0.10*rand(1,1)+1);
end
sev7_ruido=sev7_ruido';%observados con ruido
des_sev7=(0.10*rand(1,m)+1)'; %porcentaje del ruido
Cyy_1_sev7=diag(des_sev7.^-1)^2;

%Numero de capas
capas = 5;

%Resistividades
rho(1) = 30;
rho(2) = 10;
rho(3) = 10;
rho(4) = 10;
rho(5) = 60;
%incognitas
n=length(rho);
%Espesores
e(1) = 2.50;
e(2) = 6.50;
e(3) = 27.00;
e(4) = 39.0;
    
%Respuesta al modelo inicial, Sol. problema directo 
Resp=rhoa(ab,rho,e,capas,m)';

%Imprimir en columnas
fprintf(' AB/2    Ressistividades\n') 
for ii = 1:m
	rabi = ab(ii);
        fprintf('%6.2f %11.2f\n', rabi,Resp(ii)) 
end

per=0.0025;
h= rho*per;

%Calculo del Jacobiano
x=[1:length(rho)];

J=Jacobiano( rho, h, Resp, ab, e, capas, m )
figure('name','Sensibilidad J de Xmc_0')
plot(x,J,'-*'), xlabel('Capas modelo (n)'),ylabel('Datos (m)')
title('Sensibilidad J_{Xmc_{0}}')
%pause(2)
disp(' '), disp(' ')
%An??lisis de estabilidad
DET_JtJ=det(J'*J), Rango_Jacobiano=rank(J'*J), DVS_JtJ=svd(J'*J);, 
Coef=DVS_JtJ(1)/DVS_JtJ(5)
%pause(2)
%Jacobiano en Logaritmo
% disp(' '), disp(' ')
% for j=1:n
%     for i=1:m
% J_norm(i,j)=J(i,j)*(rho(j)/des_sev7(i));
%     end pause(1)
% end
% J_norm
% figure('name','Jacobiano Reescalado')
% plot(x,J_norm,'-*'),xlabel('Capas modelo (n)'),ylabel('Datos (m)')
% %pause(2)
% disp(' '), disp(' ')
% %An???lisis de estabilidad
% DET_JtJ_norm=det(J_norm'*J_norm),Rango_Jacobiano_norm=rank(J_norm'*J_norm)
% DVS_JtJ_norm=svd(J_norm'*J_norm);, Coef_norm=DVS_JtJ_norm(1)/DVS_JtJ_norm(5)
% %pause(2)
% disp(' '), disp(' ')
% 
            
%%%%%% Graficas de datos, Resp y modelo inicial %%%%%%%

  z=[e(1) (e(1)+e(2)) (e(1)+e(2)+e(3)) (e(1)+e(2)+e(3)+e(4)) (e(1)+e(2)+e(3)+e(4)+50) ];
 figure('name','Estimador Minimo Cuadrado Estandar')
 errorbar(ab,sev7,des_sev7,'xk'),set(gca,'xscale','log'),set(gca,'yscale','log'),hold on,
 loglog(ab,Resp,'r-'), hold on, grid on, xlabel('AB/2'),ylabel('\rho_{a}'),hold on
 stairs(z,rho,'k-'), legend('Dat_{obs}','Y_{X_{0}}','F_{0}')
 pause(1)
                          
 
                        
                       %%%%%%Estimadores minimo cuadrado%%%%%
                        
                        
%Estandar
Xmc_std=rho'; %solo en la primera iter entra rho, despues ya se actualiza
J_std=J;
for i=1:3;
  disp(['* * * * * * * * Iteracion: ', num2str(i),'* * * * * * * * *']);
  disp(' ')
  RespXmc=rhoa(ab,Xmc_std,e,capas,m)';
  Xmc_std= Xmc_std + inv(J_std'*J_std) * J_std' * (sev7 - RespXmc),disp(' ')
  Y_Xmc_std = rhoa(ab,Xmc_std,e,capas,m)'
  J_std= Jacobiano( Xmc_std, h, Y_Xmc_std, ab, e, capas, m );
  pause(2)
  loglog(ab,Y_Xmc_std,'-')
  hold on
  pause(2)
  stairs(z,Xmc_std,'--')
  %rho=Xmc_std'; mejor Xmc_std sea igual 
end
legend('Dat_{obs}','Y(X_{0})','\rho_{X_{0}}','Y(X_{1})','\rho_{X_{1}}','Y(X_{2})','\rho_{X_{2}}','Y(X_{3})','\rho_{X_{3}}')
%graficas finales Xmc_std
figure('name','Ajuste Final Xmc')
errorbar(ab,sev7,des_sev7,'xk'),set(gca,'xscale','log'),set(gca,'yscale','log')
hold on, grid on, xlabel('AB/2'),ylabel('\rho_{a}')
loglog(ab,Resp,ab,Y_Xmc_std), legend('Dat_{obs}','Y(X_{0})','Y(Xmc_{std})')
%  figure('name','Modelo Final Xmcs')
%  stairs(z,Xmc_std,'--')
%  figure('name','Sensibilidad J de Xmc_std')
figure('name','Residuales Xmc_std'), rang=[1:m];
resi_Xmc_std = sev7-Y_Xmc_std;, O=zeros(1,m);
plot(rang, resi_Xmc_std, rang, O),xlim([0.5 18])
 
uiwait(msgbox('Continuar con la estimaci???n ponderada?, presione OK','Fin de Xmc estandar','help'));

 
            %%%%%%%%%%%%%%%%%% Graficas Finales %%%%%%%%%%%%%%%%%%%%%%%

%Ajustes Finales
 figure('name','Ajustes Finales')
 errorbar(ab,sev7,des_sev7,'xk'),set(gca,'xscale','log'),set(gca,'yscale','log'),hold on,
 loglog(ab,Resp,'r-',ab,Y_Xmc_std,'m--'), hold on, grid on, xlabel('AB/2'),ylabel('\rho_{a}')
 legend('Dat_{obs}','Y(x)_{0}','Y(x)_{std}'),title('Xmc_{std}')

 %Modelos Finales 
 z=z*-1;
 figure('name','Modelos Finales')
 stairs(Xmc_std,z,'m'), grid on, grid minor, title('Xmc_{std}')%,xlim([-100 rho(end)]), ylim([z(end) 0])
 xlabel('\rho_{cal}'),ylabel('Profundidad'),ylim([-120 0])
 
 %Residuales Finales       
 figure('name','Residuales')
 plot(rang, resi_Xmc_std, rang, O,'k'),xlim([0 18])
 xlabel('m'),ylabel('Residual'),grid on, hold on, title('Xmc_{std}')
 xlim([1 m])%para que la linea de los residuales inicie casi en cero del eje x
 
 %Jacobianos
 figure('name','Matriz de Sensibilidad J')
 plot(x,J_std,'-*'),xlabel('Capas modelo (n)'),ylabel('Datos (m)')
 title('J_{X_{Std}}')
 