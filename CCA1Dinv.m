clc; clear; close all; disp('* * * * Inversion de datos CCA * * * *'); disp(' ')
 
% % % % % % % % LECTURA DE DATOS Y PARAMETROS DE ENTRADA % % % % % % % % %
[file, NoEst, NoReg, LonReg, Dt, W, Tras, r] = reAdfiLe('registros.dat');            
[M,nXven,f,fs,F1] =  Observados(file,NoEst,W,Dt,Tras,LonReg,NoReg);
% % % % % % % % % % % % % FRECUENCIAS A INVERTIR % % % % % % % % % % % % % 
 [finv,M2,OBS,F2] = BandaINV(fs, nXven, f, M, F1);
% % % % % % % % % % % % % % MODELADO DIRECTO % % % % % % % % % % % % % % % 
      V0 = 800;                       %m/s  
      Dv = 10;                        %m/s
   sigma = 50;                        
      Vp = V0 + Dv*exp((-finv.^2)./sigma);
     PAR = length(Vp); TPSDR = DirectoCCA(finv,r,Vp)';    
      [F3, Opcion, anss] = Fig3( finv, Vp, TPSDR, r, F1, F2);
% % % % % % % % % % % % % % % % % MODELADO INVERSO % % % % % % % % % % % % % 
     per = 0.025;                     %Perturbacion en el Jacobiano

% [Vpcal_s, Vpcal_r, Vpcal, Vpcal_rn] = ...
%     inv_process (Opcion, anss, finv, r, Vp, OBS, PAR, per, M2, TPSDR, V0, Dv, sigma, F1, F3);

switch anss
    case 'Yes'
        if Opcion == 1
                Vpcal_s = INVy(finv, r, Vp, OBS, PAR, per, M2, TPSDR, V0, Dv, sigma);
        elseif Opcion == 2
                Vpcal_r =...
                    INVtikhonov(finv, r, Vp, OBS, PAR, per, M2, TPSDR, V0, Dv, sigma);
        end                     
    case 'No'
        if Opcion == 1
            Vpcal = INVn(finv, r, Vp, OBS, PAR, per, M2, TPSDR, V0, Dv, sigma, F1, F3);
        elseif Opcion == 2
            %Hayq ue modificar esta rutina para que lea el nuevo modelo
            %inicial
            Vpcal_rn =...
                INVtikhonov_n(finv, r, Vp, OBS, PAR, per, M2, TPSDR, V0, Dv, sigma);
        end
end


%Falta terminar la funci?n Vpcal_rn 


%En lugar de hacer 4 funciones dos con si y dos con "NO" hacer solo 2, una
%con estandar y otra con tikhonov. Primero imprimir figuras de graficas,
%luego preguntar si el modelo directo es correcto, corregirlo hasta que sea
%el adecuado, escoger el esquema de inversion
