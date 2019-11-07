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
%   answer = questdlg('      Proceder con la Inversion?', ...
% 'Proceso Completado','Yes','No','No');
% % % % % % % % % % % % % % % % MODELADO INVERSO % % % % % % % % % % % % % 
     per = 0.025;                     %Perturbacion en el Jacobiano
% switch answer
%     case 'Yes'
%         Option=input('Inversi?n estandar (1)    Regularizaci?n de Tikhonov (2): ');
% %         switch Option
% %             case '1'
%         if Option == 1
%                 Vpcal_s = INVy(finv, r, Vp, OBS, PAR, per, M2, TPSDR, V0, Dv, sigma);
%         elseif Option == 2
%             %case '2'
%                 %Vpcal_r = INVtikhonov(finv, r, Vp, OBS, PAR, per, M2, TPSDR, V0, Dv, sigma);
%         else
%             
%         end
%                 
%             
%     case 'No'
%         Vpcal = INVn(finv, r, Vp, OBS, PAR, per, M2, TPSDR, V0, Dv, sigma, F1, F3);
% end


        disp('Modelo de Vp -> V0 + Dv*exp((-f^2)/sigma)')    
disp(['Modelo actual -> V0=',num2str(V0),' ' ';' ' ' 'Dv=',num2str(Dv),' ' ';' ' ' 'sigma=',num2str(sigma)])
   F5 = Fig5( finv, M2, TPSDR, r );
   F6 = Fig6;   
  Xmc = Vp';
    i = 0;
  RMS = 1; 
    I = diag(ones(1,PAR));
 alfa = 0.0001;
    Z = Jacobiano( finv, r, Vp, OBS, PAR, per, TPSDR );
    pause(2)
    
while(RMS > 0.05)
           RMS = sqrt(sum((M2 - TPSDR).^2)/length(M2));
    figure(3)
    hold on
    bar(i,RMS)  

             i = i+1;
        TPSDmc = DirectoCCA(finv,r,Xmc)';               %transpuesto solo para visualizacion
        Vp_reg = Xmc + inv(Z'*Z + (alfa*I)) * Z' * (M2 - TPSDmc);
       TPSDcal = DirectoCCA(finv,r,Vp_reg)';               %transpuesto solo para visualizacion
             Z = Jacobiano( finv, r, Vp_reg, OBS, PAR, per, TPSDcal );

    figure(2);
    hold on
           FF2 = loglog(finv,TPSDcal,'--r','LineWidth',1);
    legend('M_{Obs}','PSD_{0}',strcat('PSD_{iter:', num2str(i),'}'))
    pause(1.5)

% Determinante=det(Z'*Z);
% DVS_ZtZ=svd(Z'*Z);
% Ind=DVS_ZtZ(length(Z))/DVS_ZtZ(1);
% % Ind=DVS_ZtZ(1)/DVS_ZtZ(length(Z));
% 
% 
% if (Ind < 10^3)
%     disp(['Inestabilidad en la inversa de la matriz ZtZ segun el indice DVS=', ' ', num2str(Ind)])
%     break
% elseif(Determinante < 0)
%     disp(['Inestabilidad en la inversa de la matriz ZtZ. Determinante=', ' ' ,num2str(Determinante)])
% break
% else
%     
% end

      Residual = abs(sum(M2 - TPSDcal));
    if 0.05>Residual
        break 
    else
        delete(FF2)
    end
    legend('M_{Obs}','PSD_{0}')

    TPSDR = TPSDcal;
    disp(['Iteracion: ',num2str(i)])
    RMS
    Residual
end
    F7 = Fig7( finv, Vp_reg);








