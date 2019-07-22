function [ M ] = Observados( file, NoEst, W, Tras, nV, nXven )

% Fun_M     Pograma CCA mejorado para n archivos y ventaneado. El algoritmo                                   %
%           o flujo seguido para programarlo fue el sugerido por:          
%           Cho et al,(2004): A new method to dtermine phase velocities of 
%           Rayleigh waves from microseismis                
%                                                                          
% Creado por FCH 2009 y modificado por MAOG 2019                                                                                                                  %  
%                                                                          
% % % % % % % % % % % % Lista de Variables % % % % % % % % % % % % %
%
%       file  =    Lectura del archivo .dat que se graba en campo
%      NoEst  =    Numero de estaciones
%      NoReg  =    Numero de registos de ruido sismico
%      LonReg =    Longitud del registro de ruido en segundos
%          Dt =    Muestreo en segundos
%           W =    Ancho de la ventana en segundos
%        Tras =    Porcentaje del traslape
% 

% clc; clear all
% file=load ('registros.dat');
% NoEst=input('Numero de estaciones en el arreglo circular: ')
% NoReg=input('Numero de registros de ruido: ')
% LonReg=input('Longitud de registro en segundos: ') 
% Dt=input('Muestreo: ')
% W=input('Tama?o de la ventana en segundos: ')
% Tras=input('Traslape entre ventanas 1(0%) o 2(50%): ')
% 
% 
% 
% 
%      nV = (LonReg*NoReg)/W;      %numero de ventanas totales
%   nXven = W/Dt                   %datos por ventana             RENOMBRAR POR Nven al parecer siempre ser? igual a las observaciones por lo que se puede quitar y ahorrar una variable
%      fs = 1/Dt;                  %Frecuencia Maxima
%      ds = 1/(nXven*Dt);          %Muestreo frec
%       f = (0:ds:nXven-1)*(fs/nXven); % Frequency range

pi=4*atan(1);

%si se escoge la opcion de traslape al 50% se debe bajar el muestreo de
%0.008 a 0.004

if Tras == 2     
    nV=(nV*2)-1;
end


dteta=2*pi/NoEst;                  %delta de teta
l=1;                            %Variable que ontrola el salto de lineas en los registros
for h=1:nV                    %controla el numero de ventanas
    for jj=0:nXven-1               %controla el numero de datos por ventana 
     alpha0(jj+1)=0;            %Series de tiempo 
     alpha1(jj+1)=0;
        for k=1: NoEst            %controla el numero de estaciones
            alpha0(jj+1)=alpha0(jj+1)+file(jj+l,k)*dteta;                              %Serie de tiempo 1 con funcion de peso
            alpha1(jj+1)=alpha1(jj+1)+(file(jj+l,k)*exp(-i*1*k*2*pi/NoEst))*dteta;        %Serie de tiempo 2 con funcion de peso
        end
     alpha0n(jj+1,h)=alpha0(jj+1);   %aqui es guardar las series de tiempo obtenidas alpha0 y alpha1 ventaneadas en columnas 
     alpha1n(jj+1,h)=alpha1(jj+1);
    end
    if Tras == 2
        l=l+nXven/2;
    else
        l=l+nXven;
    end
end    
    
%Obteniendo y aplico la Ventana de Hanning a cada ventana y lo guardo en matriz
windowHanning = window(@hann,nXven).';

for ii=1:nV
    alpha0nV(:,ii)=windowHanning'.*alpha0n(:,ii);
    alpha1nV(:,ii)=windowHanning'.*alpha1n(:,ii);
end

%Se obtiene la ventana de parzen en tiempo
windowParzen=window(@parzenwin,nXven).';
%Se obtiene la norma de la ventana de parzen
FTwindowParzen=fft(windowParzen,nXven);
AbsParzen=abs(FTwindowParzen);


%Calculo de la densidad espectral de potencia G0 y G1 la forma de calculralos en ciclos FOR tambien equivale 
%utilizando indices de matrices y vectorteros eso reuce el numero de linea y el tiempo de calculo: por ejemplo G0(:,:)=fft
% 
% fs=1/Dt;
% f = (0:n-1)*(fs/n);     % Frequency range
for jj=1:nV
   G0(:,jj) = fft(alpha0nV(:,jj),nXven);                %Se calculan las FFT de cada ventana claculando primero la densidad espectral de potencia
   G1(:,jj) = fft(alpha1nV(:,jj),nXven);                % G0 y G1, 
   absG0(:,jj) =abs(G0(:,jj));                      % despues se calcula |G0| y |G1| de cada ventana
   absG1(:,jj) =abs(G1(:,jj));
   PSD_G0(:,jj)=(absG0(:,jj).^2)./W;                %Se calculan la densidad espectral de potencia |FFT|^2/V, la cual es 
   PSD_G1(:,jj)=(absG1(:,jj).^2)./W;                %el termino central ecuaci??n 16 del archivo -El metodo CCA-: Primera Reunion de Avances 
   SuavisadoPSD_G0(:,jj)=AbsParzen'.*PSD_G0(:,jj);  %Parzen
   SuavisadoPSD_G1(:,jj)=AbsParzen'.*PSD_G1(:,jj);                                         
end 

%Promedio de las densidades espectrales de potencia
SumPSD_G0=sum(PSD_G0');
SumPSD_G1=sum(PSD_G1');
PromPSD_G0=SumPSD_G0./nV;
PromPSD_G1=SumPSD_G1./nV;


%Calculamos la funcion M que es la representaci??n de los datos observados 
M = PromPSD_G0./PromPSD_G1;


%Promedio de las densidades espectrales de potencia suavisadas con Parzen
SumSuavisadoPSD_G0=sum(SuavisadoPSD_G0');
SumSuavisadoPSD_G1=sum(SuavisadoPSD_G1');
PromSuavisadoPSD_G0=SumSuavisadoPSD_G0./nV;
PromSuavisadoPSD_G1=SumSuavisadoPSD_G1./nV;
M_smooth= PromSuavisadoPSD_G0./PromSuavisadoPSD_G1; %asi se corrigio


%matrix(:)=M;
%M=matrix;
%save sal75.dat matrix -ascii


 % % % % % % Estas graficas ya estan en el programa InversionCCA.m % % % % % % %
%Graficado del cociente G0/G1 con y sin suavizado
% figure('name','Funcion M')
% 
% loglog (f,M,'b',f,M_smooth,'r');title('M[rk(w)]')
% legend('M','M_{Smooth}')
% grid on; xlabel('f'); ylabel('M')
% 
% figure('name','Funcion M')
% loglog(f,M,'b'); grid on
% 
% figure('name','Funcion M suavizada')
% loglog(f,M_smooth,'r')
% %FIN DE LAS GRAFICAs comentadas

end


