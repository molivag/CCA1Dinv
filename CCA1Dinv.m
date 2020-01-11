clc; clear; close all; disp('* * * * Inversion de datos CCA * * * *'); disp(' ')

%% % % % % % % % LECTURA DE DATOS Y PARAMETROS DE ENTRADA % % % % % % % % %
[file, NoEst, NoReg, LonReg, Dt, W, Tras, r] = reAdfiLe('registros.dat');            
[M,nXven,f,fs,F1] =  Observados(file,NoEst,W,Dt,Tras,LonReg,NoReg);

%% % % % % % % % % % % % % FRECUENCIAS A INVERTIR % % % % % % % % % % % % % 
 [finv,M2,OBS,F2] = BandaINV(fs, nXven, f, M, F1);

%% % % % % % % % % % % % % % MODELADO DIRECTO % % % % % % % % % % % % % % % 
      V0 = 800;                             %m/s  
      Dv = 10;                              %m/s
   sigma = 50;                        
      Vp = V0 + Dv*exp((-finv.^2)./sigma);  %m/s
     PAR = length(Vp); TPSDR = DirectoCCA(finv,r,Vp)';    
      [F3, Opcion] = Fig3( finv, Vp, TPSDR, r, F1, F2,V0,Dv,sigma);

%% % % % % % % % % % % % % % % % % MODELADO INVERSO % % % % % % % % % % % % 
     per = 0.025;
     VpCal= ModInv(finv, r, Vp, OBS, PAR, per, M2, TPSDR, Opcion);
     
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 