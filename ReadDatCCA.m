clc;close all; clear
% MRead_fopen = fopen('registros.dat');
MRead_dlm=dlmread('registros.dat')';

MRead_dlm2=dlmread('1-100.dat');

plot(MRead_dlm)

% Pwelch_FO=pwelch(MRead_fopen);
Pwelch_DLM=pwelch(MRead_dlm);
Pwelch_DLM2=pwelch(MRead_dlm2);

figure('name','Funcion Pwelch Registro loglog')
loglog(Pwelch_DLM); title('Registros.dat')

figure('name','Funcion PWelch 1_100 LogLog')
loglog(Pwelch_DLM2); title('1-100.dat')
