function [ Xmc, F7, F5, F6 ] = ModInv( finv, r, Vp, OBS, PAR, per, M2, TPSDR, Opcion )
     
switch Opcion
    case 1
       
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
pause(1)

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

if RMS <= 0.05
    break 
else
    delete(FF2)
end
legend('M_{Obs}','PSD_{0}')

TPSDR=TPSDcal;
% disp(['Iteracion: ',num2str(i)])


end
disp(['Iteracion: ',num2str(i)])
F7 = Fig7( finv, Xmc);
RMS=RMS*100;
disp(['RMS = ',num2str(RMS),'%'])
                    
    case 2
      
        
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


end

