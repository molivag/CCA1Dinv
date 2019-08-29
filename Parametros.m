function [file, NoEst, NoReg, LonReg, Dt, W, Tras, r] = reAdfiLe( registos )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


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

end

