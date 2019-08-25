clc; clear; close all; disp('* * * * Inversion de datos CCA * * * *'); disp(' ')

% % % % % % % % % % DATOS OBSERVADOS: PARAMETROS DE ENTRADA % % % % % % % % 
  file = load ('registros.dat'); %registros.dat antes lamado todas.dat
%NoEst=input('Numero de estaciones en el arreglo circular: ');
 NoEst = 20;
%NoReg=input('Numero de registros de ruido: ');
 NoReg = 20;
%LonReg=input('Longitud de registro (seg): ');
LonReg = 65;
%Dt=input('Muestreo: ');
    Dt = 0.004;
%W=input('Tamanio de la ventana (seg): ');
     W = 7.5;
%Tras=input('Traslape de ventanas 1(0%) o 2(50%): ');
  Tras = 1;
     r = 15;             
[M,nXven,f,fs,F1] =  Observados(file,NoEst,W,Dt,Tras,LonReg,NoReg);
% % % % % % % % % % % % % FRECUENCIAS A INVERTIR % % % % % % % % % % % % % 
 [finv,M2,OBS,F2] = BandaINV(fs, nXven, f, M, F1);
% % % % % % % % % % % % % % MODELO DIRECTO % % % % % % % % % % % % % % % %      
      V0 = 800;                       %m/s  OPTIMO 1
      Dv = 10;                        %m/s
   sigma = 50;                        %OPTIMO 0.5
      Vp = V0 + Dv*exp((-finv.^2)./sigma);
     PAR = length(Vp);          %introducir en la funcion Directo CCA y ponerlo como parametro de salida para ser leido por Jacobiano
    TPSD = DirectoCCA(finv,r,Vp)';    %transpuesto solo para visualizacion
      F3 = Fig3( finv, Vp, TPSD, r, F1, F2);
     per = 0.025;                     %Perturbacion en el Jacobiano
  answer = questdlg('      Proceder con la Inversion?', ...
'Proceso Completado','Yes','No','No');
% % % % % % % % % % % % % % INVERSI?N % % % % % % % % % % % % % % % % % % %      
switch answer
case 'Yes'
Vpcal = INVy(finv, r, Vp, OBS, PAR, per, M2, TPSD, V0, Dv, sigma);
case 'No'
Vpcal = INVn(finv, r, Vp, OBS, PAR, per, M2, TPSD, V0, Dv, sigma, F1, F3);
end
