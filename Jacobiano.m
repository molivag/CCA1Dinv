function [ Z ] = Jacobiano( Vp, OBS, PAR, per, FM )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here



        %n x %p
Z=zeros(OBS,PAR);   %Inicializacion de la matriz Jacobiano 
h = Vp*per;          %Esto es para indicar un 0.25% de VP como perturbacion

    k=1;
    for j= 1:PAR
         if j==1
            Vp1 = Vp;
            Vp2 = Vp;
            Vp1(j) = Vp(j)+h(j);
            Vp2(j) = Vp(j)+2*h(j);
            Fper_1 = DirectoCCA(f, r, Vp1, OBS)';                      %respuesta f(x+h)
         Fper_1_2h = DirectoCCA(f, r, Vp2, OBS)';                  %respuesta f(x+2h)
            for i=1:OBS
                Z(i,j)=(-3*FM(i) + 4*Fper_1(i) - Fper_1_2h(i))/(2*h(j));
            end
            
         elseif j==(PAR)
               Vp_p5 = Vp;
              Vp_2p5 = Vp;
            Vp_p5(j) = Vp(j)-h(j);%rho - perturbacion 5
           Vp_2p5(j) = (Vp(j)-(2*h(j)));%rho - 2veces perturbacion 5
              Bper_5 = DirectoCCA(f, r, Vp_p5, OBS); %respuesta f(x-h)
           Bper_5_2h = DirectoCCA(f, r, Vp_2p5, OBS); %respuesta f(x-2h)
            for i=1:OBS
                Z(i,j)=( 3*FM(i) - 4*Bper_5(i) + Bper_5_2h(i) )/ (2*h(j));
            end

         else
            k=k+1;
            Vp_pluss_centro = Vp;
            Vp_minus_centro = Vp;
         Vp_pluss_centro(k) = Vp(k)+h(k);
         Vp_minus_centro(k) = Vp(k)-h(k);
            Fper_2= DirectoCCA( f, r, Vp_pluss_centro, OBS )' ;    %FM con f(x+h)
            Bper_2= DirectoCCA( f, r, Vp_minus_centro, OBS )';
            for i=1:OBS
                Z(i,j)=(Fper_2(i) - Bper_2(i))/ (2*h(j));
            end
         end
    end

    Z;

end

