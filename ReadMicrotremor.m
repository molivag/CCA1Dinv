%Programa que concatena archivos ASCI de registros de ruido sismico con el
%equipo GEODE
%Los archivos de entrada son en ASCI y por cada registro en campo se
%generan n archivos de acuerdo al numero de canales, esto implica tener
%ventanas de registro, ahora el problema consiste en concatenralos de
%manera tal que se gener un archivo parecido al de registro continuo.
%Clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Variables de entrada

file=input('Archivo inicial (solo  nombre): ');
nfile=input('Archivo final (solo  nombre): ');
x2 =input('Numero de canales: ');
L =input('Longitud en segundos: ');
dt=input('Muestreo en segundos: ');
ini =input('Cuantas lineas tienen de encabezado los archivos: ')
%file=100; nfile=119;
%ncanales=.024;
%L=65; dt=.004; ndat=L/dt;
%ini=36

%Variables de inicio

ndat=L/dt;
y1=1;y2=ndat;
x1=1;
ncanales=x2/1000;

for i=file:nfile;
    Archivo=num2str(i);
    k=1;
    for j=.001:0.001:ncanales;
    Ext=num2str(j,'%2.3f');%con formato de impresion pareceido al que se utiliza en Fortran
    Ext = strrep(Ext, '0.', '.');
    Ruido=[Archivo,Ext];%Ruido=strcat(Archivo,Ext) Es valido a la linea anterior
    [y(:,k)] = textread(Ruido,'%f','headerlines',ini);
    k=k+1;
end;
   Registros(y1:y2,x1:x2)=y(:,:);
   y1=y2+1; y2=y2+ndat;
end 
        
%save todas.dat Registros -ASCII funciona, pero para el load me marco error

dlmwrite('1_100.dat', Registros, 'delimiter', '\t','precision', '%.6f')

