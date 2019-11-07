
clc; clear; close all; disp('* * * * Inversion de datos CCA * * * *'); disp(' ')
 
% % % % % % % % LECTURA DE DATOS Y PAR?METROS DE ENTRADA % % % % % % % % %
[file, NoEst, NoReg, LonReg, Dt, W, Tras, r] = reAdfiLe('registros.dat');            
[M,nXven,f,fs,F1] =  Observados(file,NoEst,W,Dt,Tras,LonReg,NoReg);
% % % % % % % % % % % % % FRECUENCIAS A INVERTIR % % % % % % % % % % % % % 
 [finv,M2,OBS,F2] = BandaINV(fs, nXven, f, M, F1);
% % % % % % % % % % % % % % MODELADO DIRECTO % % % % % % % % % % % % % % % 
      V0 = 800;                       %m/s  OPTIMO 1
      Dv = 10;                        %m/s
   sigma = 50;                        %OPTIMO 0.5
      Vp = V0 + Dv*exp((-finv.^2)./sigma);
     PAR = length(Vp);          
   TPSDR = DirectoCCA(finv,r,Vp)';    %transpuesto solo para visualizacion
      F3 = Fig3( finv, Vp, TPSDR, r, F1, F2);
%   answer = questdlg('      Proceder con la Inversion?', ...
% 'Proceso Completado','Yes','No','No');
% % % % % % % % % % % % % MODELADO INVERSO % % % % % % % % % % % % % % % % 
     per = 0.025;                     %Perturbacion en el Jacobiano
% switch answer
%     case 'Yes'
%         Vpcal = INVy(finv, r, Vp, OBS, PAR, per, M2, TPSDR, V0, Dv, sigma);
%     case 'No'
%         Vpcal = INVn(finv, r, Vp, OBS, PAR, per, M2, TPSDR, V0, Dv, sigma, F1, F3);
% end


 
 disp(' ')
 disp(' ')
 %for ii=1:5
 xp=Vp';
%  xp=zeros(PAR,1);
 Y=M2;
 arg=1:OBS; %vector para eje x de residuales
 Z = Jacobiano( finv, r, Vp, OBS, PAR, per, TPSDR );
 [U,Lambda,V]=svd(Z);
 Lambda=diag(Lambda) %los valores singulares labda es una matriz de nxn
 figure('name','Inicio de descomposicion en valores singulares de Z')
 k=-1;
 l=0;
 for i=1:PAR
     k=k+2;
     l=l+1;

    if l == 6 %| j == 12,    % matriz de figuras: 5 x 2
       figure('name','Soluciones por DVS'),
       k=1;
    end
    

    if(Lambda(i) < (1E-4)*Lambda(1)) % condicion de corte en lambda
        break %si lambda es menor que 10E-5, termina el proceso
    end
    
  UtY(i)=0;
 for j=1:OBS
     UtY(i)=UtY(i)+U(j,i)*Y(j);
 end

  Vtx(i)=UtY(i)/Lambda(i);
   for j=1:PAR
       xp(j)=xp(j) + Vtx(i)*V(j,i); %La solucion
   end
   
   %AQUI IRIA EL MODELO DirectoCCA
   %Quiza deba ir xp en lugar de Vp
 TPSDR_DVS_lam=DirectoCCA( finv, r, xp )';%Respuesta al modelo
 

 %Grafica del estimador DVS
 x=(1:length(Vp));
 subplot(5,2,k)
 plot(x,xp,'b',x,Vp','r'),grid on,
 title(['Solucion para la componente ',num2str(i),' de la DVS de Z'])
 xlabel('PAR'),ylabel('$\nu_p$','interpreter','latex')
 legend('Xp','\rho_{X_{0}}')

 %estos son los residuales el arg es m
 resi_XDVS = M2 - TPSDR_DVS_lam;
 subplot(5,2,k+1)
 plot(arg,resi_XDVS,'g-'),hold on
 origen = zeros(size(resi_XDVS)); %este vector genera una linea en el origen
 plot(arg,origen,'k-'),xlabel('OBS'),ylabel('Residual')
 xlim([0.75 OBS])%para que la linea de los residuales inicie casi en cero del eje x
 grid on, title(['P = ',num2str(i)]);

 i;
 normx(i)=sqrt(xp'*xp/PAR);
 normres(i)=sqrt(resi_XDVS(i)'*resi_XDVS(i)./OBS);
 %aqui considero solo las lambdas mas representativas
 lamb(i)=Lambda(i)

    if l == 6 %Este if es para generar otra figura con 
        l=1;  %las soluciones de DVS
    end
 end
 disp('Resultados DVS para todas las lambda')
 xp %Esta es la Velocidad de fase calculada
 TPSDR_DVS_lam
 lamb = lamb'
%  Zdvs =Jacobiano( finv, r, Vp, OBS, PAR, per, TPSDR_DVS_lam );

ix=1:PAR;
 %curvas de convergencia
 figure('name', 'Grafico Picard')
 title('Grafico de Picard'),grid on
 semilogy(ix,abs(UtY),'r-.x',ix,abs(Vtx),'b--x',ix,lamb,'k-x')  
 axis([0.75 5 10E-15 10E5]), legend('|u^Ty|','|v^Tx|','\lambda')
 xlabel('Parametro i')
 
 % Curva L
 figure('name','Curva L')
 loglog(normx,normres,'b--'),grid on
%  xlabel('Norma de la Solucion'), ylabel('Norma de los Residuales')
 xlabel('Norma Residuales ||Gm - d||_2'), ylabel('Norma Soluci?n ||m||_2')
 title('Curva L para las 4 Primeras Componentes SVD de Iter: '),hold on
 loglog(normx,normres,'kd','MarkerFaceColor','g'), hold on
 for i=1:length(normx),
   text(normx(i),normres(i),num2str(i))
 end


 figure('name','Estimador por DVS')
 title('Ajuste DVS para todas las \lambda')
 loglog(finv,M2,'xk');hold on, loglog(finv,TPSDR,'r-',finv,TPSDR_DVS_lam,'b--')
 grid on, xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
 y=ylabel('$ M \left[rk(w) \right]$','FontSize', 11,'Rotation',90,'interpreter','latex')
 set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0])
 legend('M[rk(w)]','TPSDR','M[rk(w)]_{DVS}')


 %aqui selecciono la P optima para calcular Xp
 disp(' ')
 disp(' ')
 P=input('Ingrese la solucion P deseada: ')
 disp(' ')
 disp('Calculando Xp hasta P')
 Xp_DVS=zeros(PAR,1);
 [U_P,Lambdas,V_P]=svd(Z);
 Lambdas=diag(Lambdas)
 for ii=1:P
       UtY_P(ii)=0.0;
     for jj=1:OBS
         UtY_P(ii)=UtY_P(ii)+U_P(jj,ii)*Y(jj);
     end
          Vtx_P(ii)=UtY_P(ii)/Lambdas(ii);
          for jj=1:PAR
            Xp_DVS(jj)=Xp_DVS(jj)+Vtx_P(ii)*V_P(jj,ii); %esa es la solucion
          end
 % Y_X_DVSP = rhoa(ab,Xp_DVS,e,capas,m)';
  Y_X_DVSP = DirectoCCA(finv,r,Xp_DVS)';
  lamb_P(ii)=Lambdas(ii)';
  end

 resi_XDVS_P = M2 - Y_X_DVSP;
 disp(' ')
 disp(['Resultados DVS hasta P= ',num2str(P)])
 disp(' ')
 Xp_DVS
 Y_X_DVSP
 lamb_P = lamb_P'

 figure('name','Estimador DVS')
 title(['Regularizacion DVS datos observados a invertir con ...P=',num2str(P)]...
      ,'FontSize', 12,'interpreter','latex')
 loglog(finv,M2,'xk'), hold on
 loglog(finv,TPSDR,'r-'), hold on, loglog(finv,Y_X_DVSP,'-.')
 grid on, xlabel('Frecuencia (Hz)','FontSize', 11,'interpreter','latex')
 y=ylabel('$ M \left[rk(w) \right]$','FontSize', 11,'Rotation',90,...
   'interpreter','latex') ;
 set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0])

% 
% 
%
%
%
%

