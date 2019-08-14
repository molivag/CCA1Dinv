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
W = 20;
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
     F1 = Fig1(f,M);

% % % % % % % % % % % % % FRECUENCIAS A INVERTIR % % % % % % % % % % % % % 

[finv, M2, OBS, F2] = BandaINV(fs, nXven, f, M, F1);

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
  per = 0.025;                                %Perturbacion en el Jacobiano
answer = questdlg('      Proceder con la Inversion?', ...
'Proceso Completado','Yes','No','No');

switch answer
case 'Yes'

Vpcal = INV(finv, r, Vp, OBS, PAR, per, M2, TPSD );
    
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

Vpcal = INV(finv, r, Vp, OBS, PAR, per, M2, TPSD );

end

