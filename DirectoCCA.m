function [ TPSD ] = DirectoCCA( f, r, Vp )

%DirectoCC  Resuelve el proble directo de la propagacion de ondas Rayleigh 
%           que se registran en un arreglo circular sin centro midiendoi 
%           la componente vertical a partir de la cual se estimara el 
%           cociente del espectro d epotencia teorico.

%           La funcion calcula el cociente de lasa funciones dfe bessel de
%           primer y segundo orden ambas d eprimera especia.
            
%           La funcion devuelve como resultado un vector de (n x 1)
%           donde n es el numnero de datos observados, OBS


% % % % % % % % % % % % Lista de Variables % % % % % % % % % % % % %

%       f  =    Ancho de banda y/o frecuencia de muestreo
%       r  =    Radio del arreglo circular
%       Vp =    Curva de velocidad de fase (Modelo Inicial)
%      OBS =    Numero de datos observados 


muestras = length(f);

%Numero de onda K
kr(muestras)=0; %inicializacion del vector de kr, que es el produucto del
           % numero de onda y el radio del arreglo circular

for i=1:muestras
kr(i) = 2*pi*r*(f(i)./Vp(i)); %se calcula el numero de onda a partir de
end                           %la velocidad de fase propuesta (Vp)


                              
%Calculo de las funciones de Bessel de primera especie, Jn(z) para cada 
%frecuencia del ancho de banda, (f).
J0T(muestras)=0;
J1T(muestras)=0;
for i=1:muestras
J0T(i) = (besselj(0,kr(i))).^2;
J1T(i) = (besselj(1,kr(i))).^2;
end

%Se calcula el cociente entre ambas funciones de Bessel
TPSD(muestras)=0;
for i=1:muestras
TPSD(i) = J0T(i)./J1T(i);  
end

end
