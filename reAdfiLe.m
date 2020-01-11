function [file, NoEst, NoReg, LonReg, Dt, W, Tras, r] = reAdfiLe( registros )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% % % % % % % % % % DATOS OBSERVADOS: PARAMETROS DE ENTRADA % % % % % % % % 
  file = load (registros); %registros.dat antes lamado todas.dat
NoEst=input('No. de est. del arreglo circular: ');
%  NoEst = 20;
disp(' ')
NoReg=input('    Numero de registros de ruido: ');
%  NoReg = 20;
disp(' ')
LonReg=input('      Longitud de registro (seg): ');
% LonReg = 65;
disp(' ')
Dt=input('                        Muestreo: ');
%     Dt = 0.004;
disp(' ')
W=input('     Tamanio de la ventana (seg): ');
%      W = 60;
disp(' ')
Tras=input('         Traslape 1(0%) o 2(50%): ');
%   Tras = 1;
%   r = 15; 
disp(' ')    
    r = input('               Radio del arreglo: ');
disp(' ')    
   

end

