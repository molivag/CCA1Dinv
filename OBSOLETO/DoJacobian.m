clc; clear

f = (1E-1:2.0:1E1);            %Rango de frecuencia

V0 = 10;   %m/s  OPRTIMO 1
Dv = 10;   %m/s
sigma = .5; %     OPTIMO 0.5
%Vp = V0 + Dv*exp((-f.^2)./sigma);
Vp = [5 2 2 2 2];
r=3;
disp('Tamanio de la matriz')
OBS = length(f)     %n
PAR = length(Vp)    %p


per=0.025;
Vp
h = Vp*per %Esto es para indicar un 0.25% de VP como perturbacion

disp('* * * * * * * * * inicio de las perturbaciones * * * * * * * * *')




       %n x %p
Z=ones(OBS,PAR)
for j=1:PAR
    j
    m=1
    for i=1:OBS
        i
        Vp_pluss_Perturbada = Vp;
        Vp_minus_Perturbada = Vp;
  
        if j==1
            Z(i,j)=10
        elseif j==(PAR)
            Z(i,j)=20
        else
            if i==m
                Vp_pluss_Perturbada(i) = Vp(i)+h(j)    %ESTE CICLO PERMITE SUMAR LA PERTURBACION 2 CON EL PARAMETRO 2
                Vp_minus_Perturbada(i) = Vp(i)-h(j)    %ESTE PERMITE RESTARLO
            end                                        %Y ESTOS CICLOS RELLENA EL RESTO DEL VECTOR SIN SUMAR
            Fper_2= FMPSD( f, r, Vp_pluss_Perturbada, OBS )     %FM con f(x+h)
            Bper_2= FMPSD( f, r, Vp_minus_Perturbada, OBS )     %FM con f(x-h)
            Z(i,j)= (Fper_2(i) - Bper_2(i))/ (2*h(j))   %Esta instruccion calcula las derivadas centradas de 2o orden desde
                                                                              % la columna 2 hasta la p-1; p= PAR 
            m=m+1;
         end%if columnas
    end %i
end %j


    
            
            
   Z         
disp('FIN')



per=0.025;
Vp
h = Vp*per %Esto es para indicar un 0.25% de VP como perturbacion

disp('* * * * * * * * * inicio de las perturbaciones * * * * * * * * *')

Vp_pluss_Perturbada = Vp 
Vp_minus_Perturbada = Vp
for k=2:(PAR-1)
    k
    m=2
    for l=2:(OBS-1)
     l
        if l==m
            Vp_pluss_Perturbada(l) = Vp(l)+h(k)    %ESTE CICLO PERMITE SUMAR LA PERTURBACION 2 CON EL PARAMETRO 2
            Vp_minus_Perturbada(l) = Vp(l)-h(k)    %ESTE PERMITE RESTARLO
        end
            Vp_pluss_Perturbada(l) = Vp(l)         %Y ESTOS CICLOS RELLENA EL RESTO DEL VECTOR SIN SUMAR
            Vp_minus_Perturbada(l) = Vp(l)         %NI RESTAR NADA, ES DECIR, CALCULA LOS COEFICIENTES DE LAS DIFERENCIAS CENTRADAS
     m=m+1                                            %DESDE j=2 hasta j= PAR-1  
%     Vp_pluss_Perturbada;
%     Vp_minus_Perturbada;
%   Fper_2= FMPSD( f, r, Vp_pluss_Perturbada, OBS )  ;  %FM con f(x+h)
%   Bper_2= FMPSD( f, r, Vp_minus_Perturbada, OBS ) ;   %FM con f(x-h)
    end
   
end

% disp('Modelo directo mas un delta de perturbacion')
% Fper_2
% disp('Modelo directo menos un delta de perturbacion')
% Bper_2

%COEFICIENTES DIFERENCIAS CENTRADAS
rho_p2=[rho(1) rho(2)+h(2) rho(3) rho(4) rho(5)];%rho + perturbacion 2
rho_p_2=[rho(1) rho(2)-h(2) rho(3) rho(4) rho(5)];%rho - perturbacion 2

%Fper_2=rhoa(ab,rho_p2,e,capas,m); %respuesta f(x+h)
%Bper_2=rhoa(ab,rho_p_2,e,capas,m);%respuesta f(x-h)


%diferencia finita centrada O(h^2)
jaco2=zeros(m,1);
for l=1:m
       
    jaco2(l)= (Fper_2(l) - Bper_2(l))/ (2*h(2));
end


J= [jaco1 jaco2 jaco3 jaco4 jaco5];


















%    ****************calculo del jacobiano*******************
