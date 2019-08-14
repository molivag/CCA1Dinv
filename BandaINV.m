function [finv, M2, OBS, F2] = BandaINV(fs, nXven, f, M, F1)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


disp('Seleccione el ancho de banda a invertir' )
   min = input('Minima: ');
   max = input('Maxima: '); 
  finv = (min:(fs/nXven):max);        
     k = find(abs(f-min) < (fs/nXven));
    M2 = M(k:length(finv)+(k-1))';
   OBS = length(M2);
    F2 = Fig2( finv, M2, F1 );



end

