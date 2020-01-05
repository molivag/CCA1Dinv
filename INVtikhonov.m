function [ Xmc, F7, F5, F6 ] = INVtikhonov(finv, r, Vp, OBS, PAR, per, M2, TPSDR, V0, Dv, sigma)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

disp(' ')
disp(' * * * * * REGULARIZACION DE TIKHONOV * * * * * ')
disp(' ')
disp('Modelo de Vp  --> V0 + Dv*exp((-f^2)/sigma)')    
disp(['Modelo actual --> V0=',num2str(V0),' ' ';' ' ' 'Dv=',num2str(Dv),' ' ';' ' ' 'sigma=',num2str(sigma)])
   F5 = Fig5( finv, M2, TPSDR, r );
   F6 = Fig6;   
  Xmc = Vp';
    i = 0;
  RMS = 1; 
    I = diag(ones(1,PAR));
%  alfa = 0.000001;
            disp(' ')
            disp('Define alfa para la regularizacion')
            disp('alfa muy pequenio --> Minimo cuadrado')
            disp('alfa grande --------> Sobre amortiguamiento')
            disp('alfa medio ---------> Regularizado')
            alfa = input('alfa = ');
            disp(' ')
    Z = Jacobiano( finv, r, Vp, OBS, PAR, per, TPSDR );
%     pause(2)
while(RMS > 0.05)
             i = i+1;
           RMS = sqrt(sum((M2 - TPSDR).^2)/length(M2));
    figure(3)
    hold on
    bar(i,RMS)
    set(gca, 'XLim', [0.5, i+.5], 'XTick', 1:1:i)

             
        TPSDmc = DirectoCCA(finv,r,Xmc)';               %transpuesto solo para visualizacion
        Vp_reg = Xmc + inv(Z'*Z + (alfa*I)) * Z' * (M2 - TPSDmc);
       TPSDcal = DirectoCCA(finv,r,Vp_reg)';               %transpuesto solo para visualizacion
             Z = Jacobiano( finv, r, Vp_reg, OBS, PAR, per, TPSDcal );

    figure(2);
    hold on
           FF2 = loglog(finv,TPSDcal,'--k','LineWidth',1);
    legend('M_{Obs}','PSD_{0}',strcat('PSD_{iter:', num2str(i),'}'))
     pause(0.5)

%     disp(['Iteracion: ',num2str(i)])
%       Residual = abs(sum(M2 - TPSDcal));
      
    if RMS<=0.05
        break 
    else
        delete(FF2)
    end
    legend('M_{Obs}','PSD_{0}')

    TPSDR = TPSDcal;
      Xmc = Vp_reg;

end
disp(['Iteracion: ',num2str(i)])
F7 = Fig7( finv, Vp_reg);
RMS=RMS*100;
disp(['RMS = ',num2str(RMS),'%'])   
end






%     I = diag(ones(1,PAR));
%     Z = Jacobiano( finv, r, Vp, OBS, PAR, per, TPSDR );
%   Xmc = Vp';
%  alfa = 0.0001;
%   RMS = 1; 
%     i = 0;
 
%  
% while(RMS>0.05)
%     i=i+1;
%     disp(['Iteracion: ', num2str(i),' XmcR']),disp(' ')
%     RMS = sqrt(sum((M2 - TPSDR).^2)/length(M2))
% 
%         TPSDmc = DirectoCCA(finv,r,Xmc)';               %transpuesto solo para visualizacion
%         Vp_reg = Xmc + inv(Z'*Z + (alfa*I)) * Z' * (M2 - TPSDmc);
%        TPSDcal = DirectoCCA(finv,r,Vp_reg)';               %transpuesto solo para visualizacion
%          Z = Jacobiano( finv, r, Vp_reg, OBS, PAR, per, TPSDcal );
% figure(10)
% hold on    
%     loglog(finv,TPSDcal,'.-'),grid on, hold on
%     %rho = Xmc_reg';
%            Xmc = Vp_reg;
%     TPSDR=TPSDcal;
% 
%            legend('M_{Obs}','PSD_{0}',strcat('PSD_{iter:', num2str(i),'}'))
% end
% 






