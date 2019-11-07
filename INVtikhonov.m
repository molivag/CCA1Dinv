function [ Xmc, F7, F5, F6 ] = INVtikhonov(finv, r, Vp, OBS, PAR, per, M2, TPSDR, V0, Dv, sigma)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


disp(' ')
disp('Modelo de Vp ---> V0 + Dv*exp((-f^2)/sigma)')    
disp(['Modelo actual ---> V0=',num2str(V0),' ' ';' ' ' 'Dv=',num2str(Dv),' ' ';' ' ' 'sigma=',num2str(sigma)])
   F5 = Fig5( finv, M2, TPSDR, r );
   F6 = Fig6;   
    I = diag(ones(1,PAR));
    Z = Jacobiano( finv, r, Vp, OBS, PAR, per, TPSDR );
  Xmc = Vp';
 alfa = 0.0001;
  RMS = 1; 
    i = 0;
pause(2)

while(RMS>0.05)
RMS = sqrt(sum((M2 - TPSDR).^2)/length(M2));
figure(3)
hold on
bar(i,RMS)  

        i=i+1;
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

TPSDR=TPSDcal;
disp(['Iteracion: ',num2str(i)])
RMS
Residual
end

F7 = Fig7( finv, Xmc);


end