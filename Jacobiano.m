function [ Z ] = Jacobiano( finv, r, Vp, OBS, PAR, per, TPSDR )

% MATRIZ JACOBIANO: Esta funciOn calcula la matriz de jacobiano para el 
%                   caso general de m datos por n parametros. Esta restrin-
%                   gido a que el numero de frecuencias sea mayor o igual 
%                   a las observaciones.
%
% % % % % % % % % % % % % Lista de Variables % % % % % % % % % % % % % % % 
%           -f 	
%           -r   	
%           -Vp 
%           -OBS 
%           -PAR
%           -per
%           -TPSD 

% Explicaci?n detallada.
%       
%           Emplea diferencias finitas centradas O(h^2) para todos los puntos 
%           centrales y para los extremos emplea diferencia finita
%           adelantada (primera columna) y atrasada (ultima columna) ambas 
%           de orden O(h^2), esta funci?n trabaja con la funcion TPSD
%           resolviendo el problema directo en cada iteraci?n.



        %m x %n
Z=zeros(OBS,PAR);   %Inicializacion de la matriz Jacobiano 
h = Vp*per;         %Esto es para indicar un 0.25% de VP como perturbacion

    k=1;
    for j= 1:PAR            %Primera columna diferencia adelantada
         if j==1
            Vp1 = Vp;
            Vp2 = Vp;
            Vp1(j) = Vp(j)+h(j);
            Vp2(j) = Vp(j)+2*h(j);
            Fper_1 = DirectoCCA(finv, r, Vp1 )';                      %respuesta f(x+h)
         Fper_1_2h = DirectoCCA(finv, r, Vp2)';                  %respuesta f(x+2h)
            for i=1:OBS
                Z(i,j)=(-3*TPSDR(i) + 4*Fper_1(i) - Fper_1_2h(i))/(2*h(j));
            end
            
         elseif j==(PAR)   %Ultima columna diferencia atrasada
               Vp_p5 = Vp;
              Vp_2p5 = Vp;
            Vp_p5(j) = Vp(j)-h(j);                  %rho - perturbacion 5
           Vp_2p5(j) = (Vp(j)-(2*h(j)));            %rho - 2veces perturbacion 5
              Bper_5 = DirectoCCA(finv, r, Vp_p5);     %respuesta f(x-h)
           Bper_5_2h = DirectoCCA(finv, r, Vp_2p5);    %respuesta f(x-2h)
            for i=1:OBS
                Z(i,j)=( 3*TPSDR(i) - 4*Bper_5(i) + Bper_5_2h(i) ) / (2*h(j));
            end

         else           %Columnas centrales diferencia centrada
            k=k+1;  
            Vp_pluss_centro = Vp;
            Vp_minus_centro = Vp;
         Vp_pluss_centro(k) = Vp(k)+h(k);
         Vp_minus_centro(k) = Vp(k)-h(k);
            Fper_2= DirectoCCA( finv, r, Vp_pluss_centro )';    %FM con f(x+h)
            Bper_2= DirectoCCA( finv, r, Vp_minus_centro )';
            for i=1:OBS
                Z(i,j)=(Fper_2(i) - Bper_2(i))/ (2*h(j));
            end
         end
    end

    Z=Z;
    
    
    % figure('name','Modelo Directo')
% loglog(f,Resp,'k'); grid on; title('Modelo Directo')
%
% figure('name','Jacobisno')
% plot(Z); title('Matriz de Sensibilidad Z')

end

