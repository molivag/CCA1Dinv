function [ Xmc, F7, F5, F6 ] = INVn(finv, r, Vp, OBS, PAR, per, M2, TPSD)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


     F5 = Fig5( finv, M2, TPSD, r);
     F6 = Fig6;   
    Xmc = Vp';
      i = 0;
    RMS = 1; 
Z = Jacobiano( finv, r, Vp, OBS, PAR, per, TPSD );
disp(' ')
while(RMS>0.05)
RMS = sqrt(sum((M2 - TPSD).^2)/length(M2));
figure(3)
hold on
bar(i,RMS)   

           i=i+1;
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

Residual = abs(sum(M2 - TPSDcal));
if 0.05>Residual
    break 
else
    delete(FF2)
end
legend('M_{Obs}','PSD_{0}')

TPSD=TPSDcal;
disp(['Iteracion: ',num2str(i)])
RMS
Residual
end

F7 = Fig7( finv, Xmc);


end

