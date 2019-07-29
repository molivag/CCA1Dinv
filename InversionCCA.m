clear all; close all
clc
% % % % % % % % % % PARAMETROS DE ENTRADA % % % % % % %  % % % %

file=load ('registros.dat'); %registros.dat antes lamado todas.dat
%NoEst=input('Numero de estaciones en el arreglo circular: ');
NoEst= 20
%NoReg=input('Numero de registros de ruido: ');
NoReg= 20
%LonReg=input('Longitud de registro (seg): ');
LonReg= 65
%Dt=input('Muestreo: ');
Dt= 0.004
%W=input('Tamanio de la ventana (seg): ');
W= 1
%Tras=input('Traslape de ventanas 1(0%) o 2(50%): ');
Tras= 1
             
     nV = (LonReg*NoReg)/W;             %numero de ventanas totales
  nXven = W/Dt;                         %datos por ventana             
     fs = 1/Dt;                         %Frecuencia Maxima
     ds = 1/(nXven*Dt);                 %Muestreo frec
     f0 = 1/W;                          %Frecuencia Fundamental
      f = (f0:ds:nXven)*(fs/nXven);      %Frequency range
  M = Observados(file, NoEst, W, Tras, nV, nXven );
M=M';
OBS = length(M)

% % % % % % % % % % % % % % MODELO DIRECTO % % % % % % % % % % % % % % %  % 
r = 15;
A =900; B=0.9;                          %Expresion que define la forma de 
Vp = A.*f.^(-B);                        %la curva de velocidad de fase Vp


%    V0 = 80;        %m/s  OPRTIMO 1
%    Dv = 10;        %m/s
% sigma = .5;        %OPTIMO 0.5
%    Vp = V0 + Dv*exp((-f.^2)./sigma);
% 
% 
%    
   
TPSD = DirectoCCA(f,r,Vp)';               %transpuesto solo para visualizacion
PAR = length(Vp)
   
% % % % % % % % % % % % % % Figuras Modelo Inicial % % % % % % %  % % % % % %

F1 = Figuras( f, Vp, M, TPSD, r );


% % % % % % % % % % % % % % MODELO INVERSION % % % % % % %  % % % % % %
%y = Jx  ----->  x = ( JtJ )^(-1)*Jt(y-F(x))

per=0.025;                             %Perturbacion en el Jacobiano
Elementos=OBS*PAR
Z = Jacobiano( f, r, Vp, OBS, PAR, per, TPSD );
F2= Figuras2( f, M, TPSD, r );

%Legend=cell(7,1);
  Xmc = Vp'
for i =1:5
 TPSDmc = DirectoCCA(f,r,Xmc)'               %transpuesto solo para visualizacion
    Xmc = Xmc + inv(Z'*Z) * Z' * ( M - TPSDmc )    
TPSDcal = DirectoCCA(f,r,Xmc)';
      Z = Jacobiano( f, r, Xmc, OBS, PAR, per, TPSDcal );
 pause(2)
loglog(f,TPSDcal,'LineWidth',0.5);
grid on
%Legend{i+2}=strcat('PSD_', num2str(i+2));
end
%legend(Legend)

disp('Estabilidad de la matriz ZtZ')
INV_ZtZ=inv(Z'*Z);
Determinante=det(Z'*Z)
Rango = rank(Z'*Z)
DVS_ZtZ=svd(Z'*Z);
Coef=DVS_ZtZ(1)/DVS_ZtZ(length(Z))
