clc; clear; close all; disp('* * * * Inversion de datos CCA * * * *')
disp(' ')


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
W = 7.5;
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
   min = input('Minima: ');
   max = input('Maxima: '); 
  finv = (min:(fs/nXven):max);        
     k = find(abs(f-min) < (fs/nXven));
    M2 = M(k:length(finv)+(k-1))';
   OBS = length(M2);
    F2 = Fig2( finv, M2, F1 );

% % % % % % % % % % % % % % MODELO DIRECTO % % % % % % % % % % % % % % %  % 
      r = 15;
%       A = 10000; 
%       B = 0.9;                          %Expresion que define la forma de 
%      Vp = A.*finv.^(-B);                        %la curva de velocidad de fase Vp
%     PAR = length(Vp);
   V0 = 800;        %m/s  OPTIMO 1
   Dv = 10;        %m/s
sigma = 50;        %OPTIMO 0.5
   Vp = V0 + Dv*exp((-finv.^2)./sigma);
   PAR = length(Vp);
   TPSD = DirectoCCA(finv,r,Vp)';               %transpuesto solo para visualizacion
     F3 = Fig3( finv, Vp, TPSD, r, F1, F2);

answer = questdlg('      Proceder con la Inversion?', ...
'Proceso Completado', ...
'Yes','No','No');
switch answer
case 'Yes'
% % %
    

% % % % % % % % % % % % % % MODELO INVERSION % % % % % % %  % % % % % %
    per = 0.025;                                %Perturbacion en el Jacobiano


F5 = Fig5( finv, M2, TPSD, r);
     F6 = Fig6;          
    Xmc = Vp';
      i = 0;
    Residual = 1; 
Z = Jacobiano( finv, r, Vp, OBS, PAR, per, TPSD );
pause(1.5)
while(Residual>0.05)
    
           i=i+1;
    INV_ZtZ = inv(Z'*Z);
     TPSDmc = DirectoCCA(finv,r,Xmc)';               %transpuesto solo para visualizacion
        Xmc = Xmc + INV_ZtZ * Z' * ( M2 - TPSDmc );
    TPSDcal = DirectoCCA(finv,r,Xmc)';
          Z = Jacobiano( finv, r, Xmc, OBS, PAR, per, TPSDcal );
figure(2);
hold on
FF2 = loglog(finv,TPSDcal,'--r','LineWidth',1);
legend('M_{Obs}','PSD_{0}',strcat('PSD_{iter:', num2str(i),'}'))
pause(1.5)

Determinante=det(Z'*Z);
DVS_ZtZ=svd(Z'*Z);
% Ind=DVS_ZtZ(length(Z))/DVS_ZtZ(1);
Ind=DVS_ZtZ(1)/DVS_ZtZ(length(Z));


if (Ind < 10^3)
    disp(['Inestabilidad en la inversa de la matriz ZtZ segun el indice DVS=', ' ', num2str(Ind)])
    break
elseif(Determinante < 0)
    disp(['Inestabilidad en la inversa de la matriz ZtZ. Determinante=', ' ' ,num2str(Determinante)])
break
else
    
end

Residual = abs(sum(M2 - TPSDcal));
if 0.05>Residual
    break 
else
    delete(FF2)
end
legend('M_{Obs}','PSD_{0}')
% RMS = sqrt(sum((M2 - Xmc).^2))
RMS = sqrt(sum((M2 - Xmc).^2)/length(M2))
Residual
figure(3)
hold on
bar(i,RMS)





end
F7 = Fig7( finv, Xmc );

    
case 'No'
opc=2;
while(opc==2)
      disp(' ')
%       disp('Defina la Vp; considere A.*finv.^(-B)')    
%       disp(['Anterior ---> A=',num2str(A),' ' ';' ' ' 'B=',num2str(B)])
      disp('Defina la Vp; considere V0 + Dv*exp((-f^2)/sigma)')    
      disp(['Anterior ---> V0=',num2str(V0),' ' ';' ' ' 'Dv=',num2str(Dv),' ' ';' ' ' 'sigma=',num2str(sigma)])



   disp(' ')
   V0 = input('V0= '); 
   Dv = input('Dv= ');                          %Expresion que define la forma de 
sigma = input('Sigma =');
   Vp = V0 + Dv*exp((-finv.^2)./sigma);
 TPSD = DirectoCCA(finv,r,Vp)';               %transpuesto solo para visualizacion
   F8 = Fig8( finv, Vp, TPSD, r, F1, F3);
%   F8=F3;
 opc = input('El modelo inciail es correcto 1(Si), 2(No): ');
 while(opc ~= 1 && opc ~= 2)
disp(' ')
disp('Error, solo se reconoce la opcion 1 y 2. Intente de nuevo') 
opc = input('El modelo inciail es correcto 1(Si), 2(No): ');
 end
 
clc
end




% % % % % % % % % % % % % % MODELO INVERSION % % % % % % %  % % % % % %
    per = 0.025;                                %Perturbacion en el Jacobiano


F5 = Fig5( finv, M2, TPSD, r);
     F6 = Fig6;          
    Xmc = Vp';
      i = 0;
    Residual = 1; 
Z = Jacobiano( finv, r, Vp, OBS, PAR, per, TPSD );
pause(1.5) %este pause es para la primera iteracion 
           %que de tiempo de verse en lo que abre la figur
while(Residual>0.05)
    
           i=i+1;
    INV_ZtZ = inv(Z'*Z);
     TPSDmc = DirectoCCA(finv,r,Xmc)';               %transpuesto solo para visualizacion
        Xmc = Xmc + INV_ZtZ * Z' * ( M2 - TPSDmc );
    TPSDcal = DirectoCCA(finv,r,Xmc)';
          Z = Jacobiano( finv, r, Xmc, OBS, PAR, per, TPSDcal );
figure(2);
hold on
FF2 = loglog(finv,TPSDcal,'--r','LineWidth',1);
legend('M_{Obs}','PSD_{0}',strcat('PSD_{iter:', num2str(i),'}'))
pause(1.5)

Determinante=det(Z'*Z);
DVS_ZtZ=svd(Z'*Z);
% Ind=DVS_ZtZ(length(Z))/DVS_ZtZ(1);
Ind=DVS_ZtZ(1)/DVS_ZtZ(length(Z));


if (Ind < 10^3)
    disp(['Inestabilidad en la inversa de la matriz ZtZ segun el indice DVS=', ' ', num2str(Ind)])
    break
elseif(Determinante < 0)
    disp(['Inestabilidad en la inversa de la matriz ZtZ. Determinante=', ' ' ,num2str(Determinante)])
break
else
    
end

Residual = abs(sum(M2 - TPSDcal));
if 0.05>Residual
    break 
else
    delete(FF2)
end
legend('M_{Obs}','PSD_{0}')
% RMS = sqrt(sum((M2 - Xmc).^2))
RMS = sqrt(sum((M2 - Xmc).^2)/length(M2))
Residual
figure(3)
hold on
bar(i,RMS)




end

F7 = Fig7( finv, Xmc );

end

