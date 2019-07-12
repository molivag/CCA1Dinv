clc; clear
close all
f = (1E-1:1.5:1E1);            %Rango de frecuencia
V0 = 80;   %m/s  OPRTIMO 1
Dv = 10;   %m/s
sigma = .5; %     OPTIMO 0.5

Vp = V0 + Dv*exp((-f.^2)./sigma);
r=3;

disp('Tamanio de la matriz')
OBS = length(f);     %n
PAR = length(Vp);    %p
        %n x %p
Z=zeros(OBS,PAR);
dim=size(Z)
x=(1:length(Vp));

per=0.075;
h = Vp*per; %Esto es para indicar un 0.25% de VP como perturbacion

Resp=DirectoCCA(f, r, Vp, OBS)';

disp('* * * * * * * * * L a  m a t r i z   d e   J a  c o b i a n o * * * * * * * * *')
%close all
    k=1;
    for j= 1:PAR
     j  ;
         if j==1
            Vp1 = Vp;
            Vp2 = Vp;
            Vp1(j) = Vp(j)+h(j) ; 
            Vp2(j) = Vp(j)+2*h(j);
            Fper_1=DirectoCCA(f, r, Vp1, OBS)';                      %respuesta f(x+h)
            Fper_1_2h=DirectoCCA(f, r, Vp2, OBS)';                  %respuesta f(x+2h)
            for i=1:OBS
                Z(i,j)=(-3*Resp(i) + 4*Fper_1(i) - Fper_1_2h(i))/(2*h(j));
            end
            
           
         elseif j==(PAR)
            Vp_p5=Vp;
            Vp_2p5=Vp;
            Vp_p5(j)=Vp(j)-h(j);%rho - perturbacion 5
            Vp_2p5(j)= (Vp(j)-(2*h(j)));%rho - 2veces perturbacion 5
            Bper_5=DirectoCCA(f, r, Vp_p5, OBS); %respuesta f(x-h)
            Bper_5_2h=DirectoCCA(f, r, Vp_2p5, OBS); %respuesta f(x-2h) 
            
            for i=1:OBS
                Z(i,j)=( 3*Resp(i) - 4*Bper_5(i) + Bper_5_2h(i) )/ (2*h(j));   
            end
            
         else
            j;
            k=k+1;
            Vp_pluss_Perturbada = Vp;
            Vp_minus_Perturbada = Vp;
            
            Vp_pluss_Perturbada(k) = Vp(k)+h(k);
            Vp_minus_Perturbada(k) = Vp(k)-h(k);
            Fper_2= DirectoCCA( f, r, Vp_pluss_Perturbada, OBS )' ;    %FM con f(x+h)
            Bper_2= DirectoCCA( f, r, Vp_minus_Perturbada, OBS )';
            for i=1:OBS
                i;
                Z(i,j)=(Fper_2(i) - Bper_2(i))/ (2*h(j)) ;
            end
         end
    end
    
    Z
Elementos=OBS*PAR
Det_ZtZ=det(Z'*Z)
%format longG

INV_ZtZ=inv(Z'*Z)
norma=norm(INV_ZtZ)
    
figure('name','Modelo Directo')
loglog(f,Resp,'k'); grid on; title('Modelo Directo')

figure('name','Jacobisno')
plot(Z); title('Matriz de Sensibilidad Z')
