function [ TPSD ] = DirectoCCA( f, r, Vp, OBS )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here




%Numero de onda K
kr(OBS)=0;
for i=1:OBS
kr(i) = 2*pi*r*(f(i)./Vp(i));
end

%computes the Bessel function of the first kind, Jn(z) 
%for each element frequency range.
J0T(OBS)=0;
J1T(OBS)=0;
for i=1:OBS
J0T(i) = (besselj(0,kr(i))).^2;
J1T(i) = (besselj(1,kr(i))).^2;
end

%Forward Model
TPSD(OBS)=0;
for i=1:OBS
TPSD(i) = J0T(i)./J1T(i);  
end
end