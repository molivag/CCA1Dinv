clc; clear; close all; disp('* * * * Inversion de datos CCA * * * *'); disp(' ')
 
% % % % % % % % LECTURA DE DATOS Y PAR?METROS DE ENTRADA % % % % % % % % %
[file, NoEst, NoReg, LonReg, Dt, W, Tras, r] = reAdfiLe('registros.dat');            
[M,nXven,f,fs,F1] =  Observados(file,NoEst,W,Dt,Tras,LonReg,NoReg);
% % % % % % % % % % % % % FRECUENCIAS A INVERTIR % % % % % % % % % % % % % 
 [finv,M2,OBS,F2] = BandaINV(fs, nXven, f, M, F1);
% % % % % % % % % % % % % % MODELADO DIRECTO % % % % % % % % % % % % % % % 
      V0 = 800;                       %m/s  OPTIMO 1
      Dv = 10;                        %m/s
   sigma = 50;                        %OPTIMO 0.5
      Vp = V0 + Dv*exp((-finv.^2)./sigma);
     PAR = length(Vp);          
   TPSDR = DirectoCCA(finv,r,Vp)';    %transpuesto solo para visualizacion
      F3 = Fig3( finv, Vp, TPSDR, r, F1, F2);
  answer = questdlg('      Proceder con la Inversion?', ...
'Proceso Completado','Yes','No','No');
% % % % % % % % % % % % % % % % % MODELADO INVERSO % % % % % % % % % % % % % 
     per = 0.025;                     %Perturbacion en el Jacobiano
switch answer
    case 'Yes'
        Option = input('Inversion estandar (1) o Regularizacion de Tikhonov (2): ');
        if Option == 1
                Vpcal_s = INVy(finv, r, Vp, OBS, PAR, per, M2, TPSDR, V0, Dv, sigma);
        elseif Option == 2
                Vpcal_r = INVtikhonov(finv, r, Vp, OBS, PAR, per, M2,...
                                      TPSDR, V0, Dv, sigma);
        end                     
    case 'No'
            Vpcal = INVn(finv, r, Vp, OBS, PAR, per, M2, TPSDR, V0, Dv, sigma, F1, F3);
end