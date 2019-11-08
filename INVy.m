function [ Xmc, F7, F5, F6 ] = INVy(finv, r, Vp, OBS, PAR, per, M2, TPSDR, V0, Dv, sigma)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


disp(' ')
disp('Modelo de Vp ---> V0 + Dv*exp((-f^2)/sigma)')    
disp(['Modelo actual ---> V0=',num2str(V0),' ' ';' ' ' 'Dv=',num2str(Dv),' ' ';' ' ' 'sigma=',num2str(sigma)])
     F5 = Fig5( finv, M2, TPSDR, r);
     F6 = Fig6;   
    Xmc = Vp';
      i = 0;
    RMS = 1; 
Z = Jacobiano( finv, r, Vp, OBS, PAR, per, TPSDR );
disp(' ')
pause(2)
while(RMS>0.05)
          i = i+1;
        RMS = sqrt(sum((M2 - TPSDR).^2)/length(M2));
figure(3)
hold on
bar(i,RMS)
set(gca, 'XLim', [0.5, i+.5], 'XTick', 1:1:i)


           
    INV_ZtZ = inv(Z'*Z);
     TPSDmc = DirectoCCA(finv,r,Xmc)';               %transpuesto solo para visualizacion
        Xmc = Xmc + INV_ZtZ * Z' * ( M2 - TPSDmc );
    TPSDcal = DirectoCCA(finv,r,Xmc)';
          Z = Jacobiano( finv, r, Xmc, OBS, PAR, per, TPSDcal );
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
disp(['Iteracion: ',num2str(i)])
% Residual = abs(sum(M2 - TPSDcal))
RMS
if RMS <= 0.05
    break 
else
    delete(FF2)
end
legend('M_{Obs}','PSD_{0}')

TPSDR=TPSDcal;

end

F7 = Fig7( finv, Xmc);
title('Velocidad de fase por regresion lineal')


end

