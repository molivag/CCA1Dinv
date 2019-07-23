%Funcion Jacobiano Final
% Este archivo demuestra que la matriz esta bien construida para el caso general 
% de mxn siempre y cuando el modelo directo no depende de la frecuencia como
% en resistividad.

clc; clear; close all

f = (1E-1:2.0:1E1);            %Rango de frecuencia
r=15;


% % % % % % % % % % % % % % % Vp Metodo 2 % % % % % % % % % % % % % % % % %
A=1000; B=0.9;         %Expresion que define la forma de la curva de
Vp = A.*f.^(-B);       %velocidad de fase Vp


disp('Tamanio de la matriz')
OBS = length(f);    %n AQU? DEBE SER length(obs) es decir de los datos
                    % y para que no tenga error, la frecuencia debe ser 
                    %mayor o igual al numero de observaciones
                    
PAR = length(Vp);   %p
        %n x %p
Z=zeros(OBS,PAR);
dim=size(Z)
Elementos=OBS*PAR


per=0.075;
h = Vp*per; %Esto es para indicar un 0.25% de VP como perturbacion

Resp=DirectoCCA(f, r, Vp)';

disp('* * * * * * * * * L a   m a t r i z   d e   J a  c o b i a n o * * * * * * * * *')
%close all
    k=1;
    for j= 1:PAR
     j;
         if j==1
            Vp1 = Vp;       %diferencia adelantada
            Vp2 = Vp;
         Vp1(j) = Vp(j)+h(j);
         Vp2(j) = Vp(j)+2*h(j);
         Fper_1 = DirectoCCA(f, r, Vp1)';              %respuesta f(x+h)
      Fper_1_2h = DirectoCCA(f, r, Vp2)';              %respuesta f(x+2h)
            for i=1:OBS
                Z(i,j)=(-3*Resp(i) + 4*Fper_1(i) - Fper_1_2h(i))/(2*h(j));
            end


         elseif j==(PAR)
            Vp_p5=Vp;           %Diferencias atrasada
            Vp_2p5=Vp;
            Vp_p5(j)=Vp(j)-h(j);%rho - perturbacion 5
            Vp_2p5(j)= (Vp(j)-(2*h(j)));%rho - 2veces perturbacion 5
            Bper_5=DirectoCCA(f, r, Vp_p5); %respuesta f(x-h)
            Bper_5_2h=DirectoCCA(f, r, Vp_2p5); %respuesta f(x-2h)

            for i=1:OBS
                 Z(i,j)=( 3*Resp(i) - 4*Bper_5(i) + Bper_5_2h(i) )/ (2*h(j));
            end

         else
            j;
            k=k+1;
            Vp_pluss_Perturbada = Vp;   %diferencias centradas
            Vp_minus_Perturbada = Vp;
            Vp_pluss_Perturbada(k) = Vp(k)+h(k);
            Vp_minus_Perturbada(k) = Vp(k)-h(k);
            Fper_2= DirectoCCA( f, r, Vp_pluss_Perturbada)' ;    %FM con f(x+h)
            Bper_2= DirectoCCA( f, r, Vp_minus_Perturbada)';
            for i=1:OBS
                i;
                Z(i,j)=(Fper_2(i) - Bper_2(i))/ (2*h(j));
            end
         end
    end

    Z
%format longG
INV_ZtZ=inv(Z'*Z)
Det_ZtZ=det(Z'*Z)
norma=norm(INV_ZtZ)


% figure('name','Modelo Directo')
% loglog(f,Resp,'k'); grid on; title('Modelo Directo')
%
% figure('name','Jacobisno')
% plot(Z); title('Matriz de Sensibilidad Z')
