clear; close all; 
clc

% % % % % % % % % % DATOS OBSERVADOS: PARAMETROS DE ENTRADA % % % % % % % % 

file=load ('registros.dat'); %registros.dat antes lamado todas.dat
%NoEst=input('Numero de estaciones en el arreglo circular: ');
NoEst= 20;
%NoReg=input('Numero de registros de ruido: ');
NoReg= 20;
%LonReg=input('Longitud de registro (seg): ');
LonReg= 65;
%Dt=input('Muestreo: ');
Dt= 0.004;
%W=input('Tamanio de la ventana (seg): ');
W = 1
%Tras=input('Traslape de ventanas 1(0%) o 2(50%): ');
Tras= 1;
             
     nV = (LonReg*NoReg)/W;             %numero de ventanas totales
  nXven = W/Dt;                         %datos por ventana             
     fs = 1/Dt;                         %Frecuencia Maxima
     f0 = 1/W;                          %Frecuencia Fundamental
      f = (f0:nXven)*(fs/nXven);        %Frequency range
     %(fs/nXven) esto es el delta f
     fn = 1/(2*Dt);
      M = Observados(file, NoEst, W, Tras, nV, nXven );
     %f2 = (f0:ds:fn);
     F1 = Fig1(f,M);

% % % % % % % % % % % % % FRECUENCIAS A INVERTIR % % % % % % % % % % % % % 

disp('Seleccione el ancho de banda a invertir' )
disp(' ')
   min = input('Frecuencia minima: ');
   max = input('Frecuencia maxima: '); 
  finv = (min:(fs/nXven):max);    
    M2 = M(min:max)';
    %M2 = M(1:length(finv))';
    OBS = length(M2);
    F2 = Fig2( finv, M2, F1 );

% % % % % % % % % % % % % % MODELO DIRECTO % % % % % % % % % % % % % % %  % 
      r = 15;
      A = 1000000; 
      B = 0.9;                          %Expresion que define la forma de 
     Vp = A.*finv.^(-B);                        %la curva de velocidad de fase Vp
    PAR = length(Vp);
%    V0 = 80;        %m/s  OPTIMO 1
%    Dv = 10;        %m/s
% sigma = .5;        %OPTIMO 0.5
%    Vp = V0 + Dv*exp((-f.^2)./sigma);
   TPSD = DirectoCCA(finv,r,Vp)';               %transpuesto solo para visualizacion
     F3 = Fig3( finv, Vp, TPSD, r, F1, F2);
uiwait(msgbox('      Iniciar Inversion?','help'));

% % % % % % % % % % % % % % MODELO INVERSION % % % % % % %  % % % % % %
    per = 0.025;                                %Perturbacion en el Jacobiano


F5 = Fig5( finv, M2, TPSD, r);
     F6 = Fig6;          
    Xmc = Vp';
      i = 0;
    res = 1; 
Z = Jacobiano( finv, r, Vp, OBS, PAR, per, TPSD );
while(res>0.05)
    
           i=i+1;
     TPSDmc = DirectoCCA(finv,r,Xmc)';               %transpuesto solo para visualizacion
        Xmc = Xmc + inv(Z'*Z) * Z' * ( M2 - TPSDmc );
    TPSDcal = DirectoCCA(finv,r,Xmc)';
          Z = Jacobiano( finv, r, Xmc, OBS, PAR, per, TPSDcal );
figure(2);
hold on
FF2 = loglog(finv,TPSDcal,'--r','LineWidth',1);
legend('M_{Obs}','PSD_{0}',strcat('PSD_{iter:', num2str(i),'}'))
pause(1.5)
res = abs(sum(M2 - TPSDcal));
if 0.05>res
    break 
else
    delete(FF2)
end
legend('M_{Obs}','PSD_{0}')
% RMS = sqrt(sum((M2 - Xmc).^2))
RMS = sqrt(sum((M2 - Xmc).^2)/length(M2))
figure(3)
hold on
bar(i,RMS)


end
% RMS

disp('Estabilidad de la matriz ZtZ')
INV_ZtZ = inv(Z'*Z);
Determinante=det(Z'*Z)
Rango = rank(Z'*Z)
DVS_ZtZ=svd(Z'*Z);
Coef=DVS_ZtZ(1)/DVS_ZtZ(length(Z))


F7 = Fig7( finv, Xmc );
