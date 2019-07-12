  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                     *** TEORÍA DE INVERSIÓN DE DATOS ***                       %   
  %                                                                                %
  %      Proyecto de clase: Solución al problema inverso del SEV                   %      
  %                                                                                %
  %  Programa que resuelve el problema directo para el Sondeo Eléctrico Vertical   %
  %  generando una curva sintetica de resistividad aparente a partir de definir    %
  %  parametros de aberturas en AB/2, numero de capas, espesores y resistividades. %
  %                                                                                %
  %  Soluciona el problema directo resolviendo la integral de Stefanescu.          %
  %                                                                                %
  %                                                                                %
  %  Marco Antonio Oliva Gutiérrez                           Geofísica Aplicada    % 
  %  Modificado de Esparza/2017                                                    %
  %                                                                                %
  %  V 5.0                                                           02/08/2017    %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear, close all
%Numero de mediciones
m = 18;
%Medciones en AB/2
ab=[ 3 5 7 10 13 16 20 25 32 40 50 65 80 100 130 160 200 250];
sev3= [16 24 29 30 32 33 28 29 24 30 33 39 43 54 60 76 116 150];
%Número de capas
capas = 5;
%Resistividades
rho(1) = 100;
rho(2) = 100;
rho(3) = 100;
rho(4) = 100;
rho(5) = 100;

%incognitas
n=length(rho);

%Espesores
e(1) = 4.6;
e(2) = 15.4;
e(3) = 12.0;
e(4) = 168.0;
%Respuesta al modelo inicial, Sol. problema directo 
Resp=rhoa(ab,rho,e,capas,m);
std=(0.010*rand(m,1)+1)
Cyy_1=diag(std.^-1)^2
for i=1:m
Resp_ruido(i)= Resp(i)*(0.010*rand(1,1)+1);
end
Resp_ruido
fprintf(' AB/2    Ressistividades\n') 
for i = 1:m
	rabi = ab(i);
        fprintf('%6.2f %11.2f\n', rabi,Resp(i)) 
end
%vector de perturbacon con misma para todas las capas 
h=0.75*ones(1,length(rho));
%Calculo del Jacobiano
J=Jacobiano( rho, h, Resp_ruido, ab, e, capas, m )
figure('name','Análisis de Sensibilidad')
x=[1:length(rho)];
plot(x,J)
xlabel('Capas modelo (n)'),ylabel('Datos (m)'), title('Sensibilidad Jacobiano')

%Análisis de estabilidad
DET_JtJ=det(J'*J);
Rango_Jacobiano=rank(J'*J)
pause(5)

disp(' ')
disp(' ')
format longG



% %Datos de campo y modelo inicial y primera iteracion con Xc
% figure('name','iteración 1')
% loglog(ab,sev3,'bx',ab,Resp,'g--'), set(gca,'xscale','log'), set(gca,'yscale','log')
% grid on, %axis([1E0 1E3 1E0 1E3])
% xlabel('AB/2')
% ylabel('Resistividad Aparente')
% title('Curva de Resistividad Aparente')

%
Xmc_std=rho' + inv(J'*J)*J'*(sev3-Resp_ruido)'
%Resp_calc_Xmc=J*Xmc
Resp_calc_Xmc=rhoa(ab,Xmc_std,e,capas,m)'
% hold on
% loglog(ab,abs(Resp_calc_Xmc),'--')



disp('Inicio de descomposiciión en valores singulares')
pause(3)
% Vector de incógnitas es: x=(1:length(rho))
% Vector de datos "verdadero" y con ruido.
Y=Resp_ruido;
Y=Y';
arg=1:m;
xp=zeros(n,1);
k=-1;
figure
for l=1:1
disp(['Iteración: ', num2str(l)])

[U,lambda,V]=svd(J);
Lambda=diag(lambda);
    
for i=1:n
   k=k+2;
   if i == 5 | i == 10,    % matriz de figuras: 6 x 2
      figure
      k=1;
   end
   if(Lambda(i) < (10^(-10))*Lambda(1)) % condicion de corte en lambda
       break %si lambda es menor que 10E-10, termina el proceso 
   end
   
UtY(i)=0.0;
for j=1:m
    UtY(i)=UtY(i)+U(j,i)*Y(j);
end

 Vtx(i)=UtY(i)/Lambda(i);
  for j=1:n
      xp(j)=xp(j)+Vtx(i)*V(j,i); %esa es la solucion
  end
  
% la y calculada con DVS se obtiene incluyendo en la funcion respuesta el
% estimador xp calculado con DVS

Ysom=rhoa(ab,xp,e,capas,m)';
%ysom=(J*xp)

subplot(5,2,k)
plot(x,xp,'b-',x,x,'r-') 
grid on

%estos son los residuales el arg es m
res=Y-Ysom;

%k+1 indica que de la iteracion 1 mas 1 =2 osea siempre las en la 
%posicion derecha de la figura
subplot(5,2,k+1)
plot(arg,res,'b-')
hold on
origen = zeros(size(res)); %este vector genera una linea en el origen
plot(arg,origen,'k-')
xlim([0.75 20])%para que la linea de los residuales inicie casi en cero del eje x
grid on
title(['p = ',num2str(i)]);

i;
normx(i)=sqrt(xp'*xp/n);
normres(i)=sqrt(res'*res/m);
lamb(i)=Lambda(i);
%aqui considero solo las lambdas mas 
%representativas
   lamb(i)=Lambda(i);
end

pause(1)
Xdvs=xp
disp(' ')
Resp_cal_DVS=Ysom
pause(5)

J=Jacobiano( Xdvs, h, Resp_cal_DVS, ab, e, capas, m );
DET_JtJ=det(J'*J);
indice= svd(J'*J);
P=indice(1)/indice(5);

disp(' ')

if DET_JtJ >= 1e-07 || P>=1000 && P<=100000
disp('El Problema es Estable debido a que: ')
DET_JtJ
P
break
else
disp('El Problema es Inestable ya que:')
DET_JtJ
P
continue
end
pause(6)


end

%curvas de convergencia
figure('name', 'convergencia'),
ix=1:i;
semilogy(ix,abs(UtY),'r-.x',ix,abs(Vtx),'b--x',ix,lamb,'k-x')
axis([0.75 5 10E-15 10E5]),grid on
legend('|u^Ty|','|v^Tx|','\lambda')
title([ 'Convergencia Iter: ',num2str(l)]) 

  % Curva L
figure('name','Curva L')
loglog(normx,normres,'b--'),grid on
xlabel('Norma de la Solución'), ylabel('Norma de los Residuales')
title(['Curva L para las 4 Primeras Componentes SVD de Iter: ',num2str(l)])
hold on
loglog(normx,normres,'kd','MarkerFaceColor','g'), hold on
for i=1:length(normx),
  text(normx(i),normres(i),num2str(i))
end


disp(['El nuevo Jacobiano filtrado con DVS en la iteración',num2str(l)])
pause(4)
J
pause(4)
disp('EL estimador minimo cuadrado recursivo: ')
figure('name','Sesibilidad Jcobiano filtrado X_{dvs}')
plot(x,J),xlabel('Capas modelo (n)'),ylabel('Datos (m)')
title('Sesibilidad Jcobiano filtrado X_{dvs}')
pause(6)

figure('name','resultados Xmc e DVS')
loglog(ab,sev3,'x',ab,Resp,'--',ab,Resp_calc_Xmc,'--'),hold on
loglog(ab,Resp_cal_DVS,'--')
grid on
legend('Dat. Obs','Mod. Inicial','Inv Xmc','Inv DVS')


 % %El estimador de minimo cuadrado estandar
%  for i= 1: 15
%  Resp=rhoa(ab,Xdvs,e,capas,m);
%  J=Jacobiano( Xdvs, h, Resp, ab, e, capas, m )
%  Xmc=Xdvs + inv(J'*J)*J'*(sev3-Resp)';
%  det(J'*J)
%  indice= svd(J'*J);
%  P=indice(1)/indice(5)
%  
%  Resp_Xmc_2=rhoa(ab,Xdvs,e,capas,m);
% 
%  pause(2)
%  hold on
%  loglog(ab,Resp_Xmc_2,'--')
%  
%  Xdvs=abs(Xmc)
%  end
% 
% x=Xmc
 
 
%El estimador de minimo cuadrado recursivo
disp(' ')
disp(' ')
disp(' Estimador mínimo cuadrado regularizado')
disp(' ')
pause(3)


std_Xdvs=0.5*randn(length(Xdvs),1)+1; %el error del modelo con dvs es decir el Xdvs
C0=diag(std_Xdvs);
XmcR0=Xdvs;
XmcR=zeros(length(rho),1);
i=1;
iter=4;
alfa=0.095;

for i=1:iter
disp(['Iteración: ', num2str(i)])
pause(2)
Resp_XmcR=rhoa(ab,Xdvs,e,capas,m);
XmcR = XmcR0 + inv(J'*Cyy_1*J + alfa*C0) * J' * Cyy_1*(sev3 - Resp_XmcR)'
Ci= (alfa*C0+J'*Cyy_1*J);

Y_cal_XmcR = rhoa(ab,XmcR,e,capas,m)'
pause(5)
J= Jacobiano( XmcR, h, Y_cal_XmcR, ab, e, capas, m );
hold on
plot(ab,Y_cal_XmcR,'--'),grid on

XmcR0=XmcR;
Xdvs=XmcR;
C0=Ci;

end



figure('name','Residuales')
res_iter_1=Resp-Resp_calc_Xmc';
res_iter_2=Resp-Resp_cal_DVS';
origen=zeros(1,length(ab));
subplot(2,1,1)
plot(ab,origen,'k',ab,res_iter_1), xlabel('Datos'), ylabel('Residual')
title(['Iteración: ',num2str(l),' con X_{mcs}'])
subplot(2,1,2)
plot(ab,origen,'k',ab,res_iter_2), xlabel('Datos'), ylabel('Residual')
title(['Iteración: ',num2str(l),' con X_{DVS}'])




