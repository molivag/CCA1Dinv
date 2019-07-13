function [ J ] = Jacobiano( rho, h, Resp, ab, e, capas, m )
%MATRIZ JACOBIANO: Esta funci�n calcula la matriz jacobiano para un modelo 
%                  de 5 capas y 20 mediciones, no esta dise�ado para el
%                  caso general pero puede modificarse empleando ciclos for.
%
%           - rho, vector con las resistividades de las capas del modelo incluyendo
%                  el semiespacio, de (n x 1).
%           - h,   es un vector de (1 x n) que incluye el valor de perturbaci�n 
%                  para cada capa.
%           - rhoa, Es la funci�n respuesta de (n x 1) que resuelve el problema
%                   directo, son los valores de Rhoa, es decir, la curva 
%                   sint�tica de campo.
%           -ab2, Es un vector de (1 x n) que incluye las mediciones en AB/2
%           -e, Vector de (1 x n) con los espesores de las capas del modelo.
%           -capas, Numero de capas, para esta funcion, restringido a 5.
%           -m, N�mero de mediciones para este caso 20.

%   Detailed explanation goes here
%       
%           Emplea diferencias finitas centradas O(h^2) para todos los puntos 
%           centrales y para los extremos ocupa diferencia finita adelantada 
%           y atrasada ambas de orden O(h^2), esta funci�n trabaja
%           con ayuda de la funci�n respuesta, que es la funci�n que
%           resuelve el problema directo.

 
         
%  for i=1:length(rho)
%      h(i)= rho(i)*per;
%  end
%capa 1
rho_p1=[rho(1)+h(1) rho(2) rho(3) rho(4) rho(5)]; %rho + perturbacion 1
rho_2p1=[rho(1)+(2*h(1)) rho(2) rho(3) rho(4) rho(5)];%rho + 2perturbacion 1
Fper_1=rhoa(ab,rho_p1,e,capas,m); %respuesta f(x+h)
Fper_1_2h=rhoa(ab,rho_2p1,e,capas,m); %respuesta f(x+2h)

%capa 2
rho_p2=[rho(1) rho(2)+h(2) rho(3) rho(4) rho(5)];%rho + perturbacion 2
rho_p_2=[rho(1) rho(2)-h(2) rho(3) rho(4) rho(5)];%rho - perturbacion 2
Fper_2=rhoa(ab,rho_p2,e,capas,m); %respuesta f(x+h)
Bper_2=rhoa(ab,rho_p_2,e,capas,m);%respuesta f(x-h)


%capa 3
rho_p3=[rho(1) rho(2) rho(3)+h(3) rho(4) rho(5)];%rho capa 3 + perturbacion 3
rho_p_3=[rho(1) rho(2) rho(3)-h(3) rho(4) rho(5)];%rho capa 3 - perturbacion 3
Fper_3=rhoa(ab,rho_p3,e,capas,m); %respuesta f(x+h)
Bper_3=rhoa(ab,rho_p_3,e,capas,m);%respuesta f(x-h)


%capa 4
rho_p4=[rho(1) rho(2) rho(3) rho(4)+h(4) rho(5)];%rho capa 4 + perturbacion 4
rho_p_4=[rho(1) rho(2) rho(3) rho(4)-h(4) rho(5)];%rho capa 4 - perturbacion 4
Fper_4=rhoa(ab,rho_p4,e,capas,m); %respuesta f(x+h)
Bper_4=rhoa(ab,rho_p_4,e,capas,m);%respuesta f(x-h)


%semiespacio
rho_p5=[rho(1) rho(2) rho(3) rho(4) rho(5)-h(5)];%rho - perturbacion 5
rho_2p5=[rho(1) rho(2) rho(3) rho(4) rho(5)-(2*h(5))];%rho - 2veces perturbacion 5
Bper_5=rhoa(ab,rho_p5,e,capas,m); %respuesta f(x-h)
Bper_5_2h=rhoa(ab,rho_2p5,e,capas,m); %respuesta f(x-2h) 



%    ****************calculo del jacobiano*******************
jaco1=zeros(m,1);

%diferencia finita adelantada O(h^2)
for i=1:m%length(rho)
    
        jaco1(i)= (-3*Resp(i) + 4*Fper_1(i) - Fper_1_2h(i))/(2*h(1));
        %jaco1(i)= (Fper_1(i) - Resp(i))/h(1);
end
%diferencia finita centrada O(h^2)
jaco2=zeros(m,1);
for i=1:m
       
    jaco2(i)= (Fper_2(i) - Bper_2(i))/ (2*h(2));
end

%diferencia finita centrada O(h^2)
jaco3=zeros(m,1);
for i=1:m

    jaco3(i)= (Fper_3(i) - Bper_3(i))/ (2*h(3)); 
end

%diferencia finita centrada O(h^2)
jaco4=zeros(m,1);
for i=1:m
    
    jaco4(i)= (Fper_4(i) - Bper_4(i))/ (2*h(4));   
end

%diferencia finita atrasada O(h^2)
jaco5=zeros(m,1);
for i=1:m
    
    jaco5(i)= (3*Resp(i) - 4*Bper_5(i) + Bper_5_2h(i))/ (2*h(5));   
    %jaco1(i)= (Resp(i) - Bper_5(i))/h(5);

end

format longG
J= [jaco1 jaco2 jaco3 jaco4 jaco5];
end
