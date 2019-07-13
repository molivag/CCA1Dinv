% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Pograma CCA mejorado para n archivos y ventaneado                                   %
%   El alalpha0ritmo o flujo seguido para programarlo fue seguir el sugerido por:       %
%   Ikuo Cho, Taku Tada en el articulo: A new method to dtermine phase velocities of    % 
%                                       Rayleigh waves from microseismis                %
%                                                                                       %  
%                                                                                       %
%                                                                                       %
% Creado por FCH 2009 y modificado por MAOG 2019                                        %  
%                                                                                       %  
%                                                                                       %   
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clc
clear ; close all
% 
s=load ('registros.dat');
NE=input('Numero de estaciones en el arreglo circular: ')
NReg=input('Numero de registros de ruido: ')
L=input('Longitud de registro en segundos: ') 
dt=input('Muestreo: ')
V=input('Tama??o de la ventana en segundos: ')
Tras=input('Traslape entre ventanas 1(0%) o 2(50%): ')
pi=3.141592654;
nV=(L*NReg)/V; %numero de ventanas totales
n=V/dt; %numero de datos por ventana

if Tras == 2;
    nV=(nV*2)-1;
end;


dteta=2*pi/NE;                  %delta de teta
l=1;                            %Variable que ontrola el salto de lineas en los registros
for h=1:nV;                     %controla el numero de ventanas
    for jj=0:n-1;               %controla el numero de datos por ventana 
     alpha0(jj+1)=0;            %Series de tiempo 
     alpha1(jj+1)=0;
        for k=1: NE;            %controla el numero de estaciones
            alpha0(jj+1)=alpha0(jj+1)+s(jj+l,k)*dteta;                              %Serie de tiempo 1 con funcion de peso
            alpha1(jj+1)=alpha1(jj+1)+(s(jj+l,k)*exp(-i*1*k*2*pi/NE))*dteta;        %Serie de tiempo 2 con funcion de peso
        end;
     alpha0n(jj+1,h)=alpha0(jj+1);   %aqui es guardar las series de tiempo obtenidas alpha0 y alpha1 ventaneadas en columnas 
     alpha1n(jj+1,h)=alpha1(jj+1);
    end;
    if Tras == 2;
        l=l+n/2;
    else;
        l=l+n;
   end;
end;    
    
%Obteniendo y aplico la Ventana de Hanning a cada ventana y lo guardo en matriz
windowHanning = window(@hann,n).';

for ii=1:nV
    alpha0nV(:,ii)=windowHanning'.*alpha0n(:,ii);
    alpha1nV(:,ii)=windowHanning'.*alpha1n(:,ii);
end;

%Se obtiene la ventana de parzen en tiempo
windowParzen=window(@parzenwin,n).';
%Se obtiene la norma de la ventana de parzen
FTwindowParzen=fft(windowParzen,n);
AbsParzen=abs(FTwindowParzen);


%Calculo de la densidad espectral de potencia G0 y G1 la forma de calculralos en ciclos FOR tambien equivale 
%utilizando indices de matrices y vectorteros eso reuce el numero de linea y el tiempo de calculo: por ejemplo G0(:,:)=fft

fs=1/dt;
f = (0:n-1)*(fs/n);     % Frequency range
for jj=1:nV;
   G0(:,jj) = fft(alpha0nV(:,jj),n);                %Se calculan las FFT de cada ventana claculando primero la densidad espectral de potencia
   G1(:,jj) = fft(alpha1nV(:,jj),n);                % G0 y G1, 
   absG0(:,jj) =abs(G0(:,jj));                      % despues se calcula |G0| y |G1| de cada ventana
   absG1(:,jj) =abs(G1(:,jj));
   PSD_G0(:,jj)=(absG0(:,jj).^2)./V;                %Se calculan la densidad espectral de potencia |FFT|^2/V, la cual es 
   PSD_G1(:,jj)=(absG1(:,jj).^2)./V;                %el termino central ecuaci??n 16 del archivo -El metodo CCA-: Primera Reunion de Avances 
   SuavisadoPSD_G0(:,jj)=AbsParzen'.*PSD_G0(:,jj);  %Parzen
   SuavisadoPSD_G1(:,jj)=AbsParzen'.*PSD_G1(:,jj);                                         
end; 

%Promedio de las densidades espectrales de potencia
SumPSD_G0=sum(PSD_G0');
SumPSD_G1=sum(PSD_G1');
PromPSD_G0=SumPSD_G0./nV;
PromPSD_G1=SumPSD_G1./nV;


%Calculamos la funci??n M que es la representaci??n de los datos observados 
M= PromPSD_G0./PromPSD_G1;


%Promedio de las densidades espectrales de potencia suavisadas con Parzen
SumSuavisadoPSD_G0=sum(SuavisadoPSD_G0');
SumSuavisadoPSD_G1=sum(SuavisadoPSD_G1');
PromSuavisadoPSD_G0=SumSuavisadoPSD_G0./nV;
PromSuavisadoPSD_G1=SumSuavisadoPSD_G1./nV;
M_smooth= PromSuavisadoPSD_G0./PromPSD_G1;


%Graficado del cociente G0/G1 con y sin suavizado
figure('name','Funcion M')
loglog (f,M,'b',f,M_smooth,'r');title('M[rk(w)]')
legend('M','M_{Smooth}')
grid on; xlabel('f'); ylabel('M')

figure('name','Funcion M')
loglog(f,M,'b'); grid on

figure('name','Funcion M suavizada')
loglog(f,M_smooth,'r')

%matrix(:,1)=f;
matrix(:)=M;
Fun_M=matrix;

%save sal75.dat matrix -ascii