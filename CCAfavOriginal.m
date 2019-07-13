%Pograma CCA mejorado para n archivos y ventaneado
%El algoritmo o flujo seguido para programarlo fue seguir el sugerido por:
%Ikuo Cho, Taku Tada en el articulo: a new method to dtermine phase
%velocities of Rayleigh waves from microseismis
%clear all
cd('K:\Metodo CCA\Datos Tecoman 2008\ruido_circulo_GEODE')
s=load ('registros.dat');
NE=input('Numero de estaciones en el arreglo circular: ')
NReg=input('Numero de registros de ruido: ')
L=input('Longitud de registro en segundos: ') 
dt=input('Muestreo: ')
V=input('Tama?o de la ventana en segundos: ')
Tras=input('Traslape entre ventanas 1(0%) o 2(50%): ')
pi=3.141592654;
nV=(L*NReg)/V; %numero de ventanas totales
n=V/dt;        %numero de datos por ventana
if Tras == 2;
    nV=(nV*2)-1;
end;


dteta=2*pi/NE; %delta de teta
l=1; %Variable que ontrola el salto de lineas en los registros
for h=1:nV; %controla el numero de ventanas
    for jj=0:n-1;%controla el numero de datos por ventana 
        Go(jj+1)=0; G1(jj+1)=0;
        for k=1: NE;%controla el numero de estaciones
            Go(jj+1)=Go(jj+1)+s(jj+l,k)*dteta;
            G1(jj+1)=G1(jj+1)+(s(jj+l,k)*exp(-i*1*k*2*pi/NE))*dteta;
        end;
        Gon(jj+1,h)=Go(jj+1);%aqui es guardar las series obtenidas Go y G1 ventaneadas en columnas 
        G1n(jj+1,h)=G1(jj+1);
   end;
       if Tras == 2;
            l=l+n/2;
        else;
            l=l+n;
        end;
end;    
    
%Obteniendo y aplico la Ventana de Hanning a cada ventana y lo guardo en
%matriz
windowHanning = window(@hann,n).';

for ii=1:nV;
    GonV(:,ii)=windowHanning'.*Gon(:,ii);
    G1nV(:,ii)=windowHanning'.*G1n(:,ii);
end;

%Se obtiene la ventana de parzen en tiempo
windowParzen=window(@parzenwin,n).';
%Se obtiene la norma de la ventana de parzen
FTwindowParzen=fft(windowParzen,n);
AbsParzen=abs(FTwindowParzen)


%Calculo de la densidad espectral de energia para mabas series G0 y G1
%la forma de calulralos en ciclos FOR tambien equivale utilizando indices
%de matrices y vectorteros eso reuce el numero de linea y el tiempo de
%calculo: por ejemplo FFTGo(:,:)=fft
fs=1/dt;
f = (0:n-1)*(fs/n);     % Frequency range
for jj=1:nV;
    FFTGo(:,jj) = fft(GonV(:,jj),n);           %Se calculan las FFT de cada ventana
    FFTG1(:,jj) = fft(G1nV(:,jj),n);          
    powerGo(:,jj) =abs(FFTGo(:,jj));            %Se calculan las |FFT| de cada ventana
    powerG1(:,jj) =abs(FFTG1(:,jj));
    PSDGo(:,jj)=(powerGo(:,jj).^2)./V;              %Se calculan la denisas espectral de potencia |FFT|^2/V por ventana
    PSDG1(:,jj)=(powerG1(:,jj).^2)./V; 
    SuavisadoPSDGo(:,jj)=AbsParzen'.*PSDGo(:,jj) ; %Parzen
    SuavisadoPSDG1(:,jj)=AbsParzen'.*PSDG1(:,jj);
    %SmootGo(:,jj) = smooth(f,PSDGo(:,jj),0.03,'rloess'); %Suavisado
    %SmootG1(:,jj) = smooth(f,PSDG1(:,jj),0.03,'rloess');
%    SumGo=SumG1+SmootGo(:,jj);                                               
%    SumG1=SumG1+SmootG1(:,jj);                                               
end; 
%Promedio de las densidades espectrales de potencia
SumPSDGo=sum(PSDGo');
SumPSDG1=sum(PSDG1');
PromPSDGo=SumPSDGo./nV;
PromPSDG1=SumPSDG1./nV;

%Promedio de las densidades espectrales suavisadas con Parzen
SumSuavisadoPSDGo=sum(SuavisadoPSDGo');
SumSuavisadoPSDG1=sum(SuavisadoPSDG1');
PromSuavisadoPSDGo=SumSuavisadoPSDGo./nV;
PromSuavisadoPSDG1=SumSuavisadoPSDG1./nV;

%Promedio de las densidades espectrales de potencia suavisadas
%SumGo=sum(SmootGo')
%SumG1=sum(SmootG1')
%PromGo=SumGo./nV;
%PromG1=SumG1./nV;

%Graficado del cociente Go/G1
loglog (f,PromPSDGo./PromPSDG1,f,PromSuavisadoPSDGo./PromSuavisadoPSDG1,'r'); grid;title('M[rk(w)]');xlabel('f'); ylabel('M')
figure(1)
loglog(f,PromPSDGo./PromPSDG1)
figure(2)
loglog(f,PromSuavisadoPSDGo./PromSuavisadoPSDG1)
matrix(:,1)=f;
matrix(:,2)=PromPSDGo./PromPSDG1;
save sal75.dat matrix -ascii