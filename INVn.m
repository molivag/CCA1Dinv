function [ Xmc, F7, F5, F6 ] = INVn(finv, r, Vp, OBS, PAR, per, M2, TPSDR, V0, Dv, sigma, F1, F3)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

          opc=2;
while(opc==2)
      disp(' ')
%       disp('Defina la Vp; considere A.*finv.^(-B)')    
%       disp(['Anterior ---> A=',num2str(A),' ' ';' ' ' 'B=',num2str(B)])
%         disp(' ')
%          A = input('A = '); 
%          B = input('B = ');                          %Expresion que define la forma de 
%         Vp = A.*finv.^(-B); 

      disp('Defina la Vp; considere V0 + Dv*exp((-f^2)/sigma)')    
      disp(['Anterior ---> V0=',num2str(V0),' ' ';' ' ' 'Dv=',num2str(Dv),' ' ';' ' ' 'sigma=',num2str(sigma)])
   disp(' ')
   V0 = input('   V0 = '); 
   Dv = input('   Dv = ');                          %Expresion que define la forma de 
sigma = input('Sigma = ');
   Vp = V0 + Dv*exp((-finv.^2)./sigma);
 TPSDR = DirectoCCA(finv,r,Vp)';               %transpuesto solo para visualizacion
   F8 = Fig8( finv, Vp, TPSDR, r, F1, F3);
%   F8=F3;
 opc = input('El modelo inciail es correcto 1(Si), 2(No): ');
 while(opc ~= 1 && opc ~= 2)
disp(' ')
disp('Error, solo se reconoce la opcion 1 y 2. Intente de nuevo') 
opc = input('El modelo inciail es correcto 1(Si), 2(No): ');
 end
 
 if opc == 1
     continue
 else
clc
 end
end

     F5 = Fig5( finv, M2, TPSDR, r);
     F6 = Fig6;   
    Xmc = Vp';
      i = 0;
    RMS = 1; 
Z = Jacobiano( finv, r, Vp, OBS, PAR, per, TPSDR );
disp(' ')
pause(2)
while(RMS>0.05)
RMS = sqrt(sum((M2 - TPSDR).^2)/length(M2));
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

TPSDR=TPSDcal;
disp(['Iteracion: ',num2str(i)])
RMS
Residual
end

F7 = Fig7( finv, Xmc);


end

