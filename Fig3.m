function [ F3, Opcion] = Fig3( finv, Vp, TPSDR, r, F1, F2,V0,Dv,sigma)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

F3 = subplot(3,7,[12:14 19:21]);
semilogx(finv,Vp,'+k','LineWidth',1);grid on;
title('Velocidad de fase $\nu_p$ inicial','FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
y=ylabel('$\nu_p$','FontSize', 13,'FontWeight','bold','Rotation',0,'interpreter','latex');
set(gca,'YAxisLocation','right');
set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0])

% linkaxes([F1,F2,F3],'x')


F4 = subplot(3,7,15:18);
loglog(finv,TPSDR','color',[0.06,0.7,0.06]','LineWidth',0.5)
% title(['Modelo directo del cociente del espectro de densidad de potencia con radio =' ' ',num2str(r), ' ' ' m'], 'FontSize', 12,'interpreter','latex')
title(['Modelo directo de la funcion M con radio =' ' ',num2str(r), ' ' ' m'], 'FontSize', 12,'interpreter','latex')
xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex');grid on
y=ylabel('$\displaystyle \left[\frac{J_{_0}}{J_{_1} }\right]^2 $','FontSize', 11,'Rotation',0,'interpreter','latex');
set(y, 'Units', 'Normalized', 'Position', [-0.11, 0.4, 0])
%set(y, 'Units', 'Normalized', 'Position', [-0.08, 0.3, 0])
% set(gca,'YAxisLocation','right');
linkaxes([F1,F2,F4],'xy')


anss = questdlg('      Corregir el Modelo Directo? ', ...
'Proceso Completado','Yes','No','No');


switch anss
    case 'Yes'

  opc=2;
          disp(' ')
          disp('- - - - - - - - - Modelo Inicial - - - - - - - - ')
          disp(' ')
          disp('Defina la Vp; considere V0 + Dv*exp((-f^2)/sigma)')    
    disp(['Modelo actual ---> V0=',num2str(V0),' ' ';' ' ' 'Dv=',num2str(Dv),' ' ';' ' ' 'Sigma=',num2str(sigma)])
while(opc==2)
%       disp(' ')
%       disp('Defina la Vp; considere V0 + Dv*exp((-f^2)/sigma)')    
%        disp(['Modelo actual ---> V0=',num2str(V0),' ' ';' ' ' 'Dv=',num2str(Dv),' ' ';' ' ' 'Sigma=',num2str(sigma)])
   disp(' ')
   V0 = input('   V0 = '); 
   Dv = input('   Dv = ');                          %Expresion que define la forma de 
sigma = input('Sigma = ');
   Vp = V0 + Dv*exp((-finv.^2)./sigma);
 TPSDR = DirectoCCA(finv,r,Vp)';               %transpuesto solo para visualizacion
   F8 = Fig8( finv, Vp, TPSDR, r, F1);
 opc = input('El modelo incial es correcto 1(Si), 2(No): ');
 while(opc ~= 1 && opc ~= 2)
disp(' ')
disp('Error, solo se reconoce la opcion 1 y 2. Intente de nuevo') 
opc = input('El modelo incial es correcto 1(Si), 2(No): ');
 end
 
 if opc == 1
     continue
 else

 end
end
   
    case 'No'
       
end

disp(' ')
Opcion = input([' + + + + + + + Tipo de Inversion + + + + + + + +'...
'\n +                                             + ',...
'\n + - - - 1.- Minimo Cuadrado Estandar: - - - - + ',...
'\n + - - - 2.- Regularizacion de Tikhonov: - - - + ',...
'\n +                                             + ',...
'\n + + + + + + + + + + + + + + + + + + + + + + + + ',...
'\n Opcion = ']);
if Opcion == 1
    disp ('Inversion por Minimo Cuadrado Estandar')
elseif Opcion == 2
    disp('Inversion por Regularizacion de Tikhonov')
end




end

