              %%%%%%%%% Inicio DVS %%%%%%%%%

 disp(' ')
 disp(' ')
 figure
 disp('Inicio de descomposicii???n en valores singulares de J')
 pause(0.5)
 k=-1;
 %for ii=1:5
 %xp=rho';
 xp=zeros(n,1);
 Y=sev7;
 arg=1:m; %vector para eje x de residuales
 [U,Lambda,V]=svd(J);
 Lambda=diag(Lambda);  
 for i=1:n
    k=k+2;
    
    if i == 6 | i == 12,    % matriz de figuras: 5 x 2
       figure('name','Soluciones por DVS'),
       k=1;
    end
    
    if(Lambda(i) < (10E-4)*Lambda(1)) % condicion de corte en lambda
        break %si lambda es menor que 10E-5, termina el proceso 
    end
 UtY(i)=0.0;
 for j=1:m
     UtY(i)=UtY(i)+U(j,i)*Y(j);
 end
 
  Vtx(i)=UtY(i)/Lambda(i);
   for j=1:n
       xp(j)=xp(j) + Vtx(i)*V(j,i); %esa es la solucion
   end
 
 Y_DVS_lam=rhoa(ab,xp,e,capas,m)';%La solucion 
 
 %Grafica del estimador DVS 
 subplot(5,2,k)
 plot(x,xp,'b--',x,rho','r-'),grid on,
 title(['Solucion para la componente ',num2str(i),' de la DVS de J'])
 xlabel('n'),ylabel('\rho')
 legend('Xp','\rho_{X_{0}}')
 
 %estos son los residuales el arg es m
 resi_XDVS = sev7 - Y_DVS_lam;
 subplot(5,2,k+1)
 plot(arg,resi_XDVS,'g-'),hold on
 origen = zeros(size(resi_XDVS)); %este vector genera una linea en el origen
 plot(arg,origen,'k-'),xlabel('m'),ylabel('Residual')
 xlim([0.75 m])%para que la linea de los residuales inicie casi en cero del eje x
 grid on
 title(['P = ',num2str(i)]);
 
 i;
 normx(i)=sqrt(xp'*xp/n);
 normres(i)=sqrt(resi_XDVS'*resi_XDVS/m);
 %aqui considero solo las lambdas mas representativas
 lamb(i)=Lambda(i);
 
 end
 disp('Resultados DVS para todas las lambda')
 xp
 Y_DVS_lam
 lamb = lamb'
 Jdvs= Jacobiano( xp, h, Y_DVS_lam, ab, e, capas, m );

ix=1:i;
 %curvas de convergencia
 figure('name', 'convergencia')
 semilogy(ix,abs(UtY),'r-.x',ix,abs(Vtx),'b--x',ix,lamb,'k-x'),title('Convergencia') 
 axis([0.75 5 10E-15 10E5]),grid on
 legend('|u^Ty|','|v^Tx|','\lambda')
 
 % Curva L
 figure('name','Curva L')
 loglog(normx,normres,'b--'),grid on
 xlabel('Norma de la Soluci???n'), ylabel('Norma de los Residuales')
 title('Curva L para las 4 Primeras Componentes SVD de Iter: '),hold on
 loglog(normx,normres,'kd','MarkerFaceColor','g'), hold on
 for i=1:length(normx),
   text(normx(i),normres(i),num2str(i))
 end
 
 figure('name','Estimador por DVS')
 errorbar(ab,sev7,des_sev7,'xk'),set(gca,'xscale','log'),set(gca,'yscale','log'),hold on,
 loglog(ab,Resp,'r-',ab,Y_DVS_lam,'b--'),grid on, xlabel('AB/2'),ylabel('\rho_{a}')
 hold on, stairs(z,xp,'b'), legend('Dat_{obs}','Y(x)_{0}','Y(x)_{dvs}','\rho_{x_{dvs}}')
 title('Ajuste DVS para todas las \lambda')

 %aqui selecciono la P optima para calcular Xp
 disp(' ')
 disp(' ')
 P=input('Ingrese la soluci???n P deseada: ')
 disp(' ')
 disp('Calculando Xp hasta P')
 disp('...')
 pause(1)
 disp('.....')
 pause(1)
 disp('.......')
 pause(1)
 disp('.........')
 disp(' ')
 
 Xp_DVS=zeros(n,1);
 [U_P,Lambdas,V_P]=svd(J);
 Lambdas=diag(Lambdas)
 
 for ii=1:P
       UtY_P(ii)=0.0;
     for jj=1:m
         UtY_P(ii)=UtY_P(ii)+U_P(jj,ii)*Y(jj);
         end
          Vtx_P(ii)=UtY_P(ii)/Lambdas(ii);
     for jj=1:n
          Xp_DVS(jj)=Xp_DVS(jj)+Vtx_P(ii)*V_P(jj,ii); %esa es la solucion
    end
 
 Y_X_DVSP = rhoa(ab,Xp_DVS,e,capas,m)';
 
 lamb_P(ii)=Lambdas(ii)';
 end
 resi_XDVS_P = sev7 - Y_X_DVSP;
 disp(' ')
 disp(['Resultados DVS hasta P= ',num2str(P)])
 disp(' ')
 Xp_DVS
 Y_X_DVSP
 lamb_P = lamb_P'
 
 figure('name','Estimador por DVS')
 errorbar(ab,sev7,des_sev7,'xk'),set(gca,'xscale','log'),set(gca,'yscale','log'),hold on,
 loglog(ab,Resp,'r-'), hold on, grid on, xlabel('AB/2'),ylabel('\rho_{a}'),hold on
 stairs(z,rho,'k-'), legend('Dat_{obs}','Y(x)_{0}','\rho_{x_{0}}'), pause(1)
 hold on, loglog(ab,Y_X_DVSP,'-.'), hold on
 stairs(z,Xp_DVS,'--'), legend('Dat_{obs}','Y(x)_{0}','\rho_{x_{0}}','Y(x)_{dvs}','\rho_{x_{dvs}}')
 title(['Ajuste DVS con P=',num2str(P)])

  
uiwait(msgbox('Presione Ok para imprimir las graficas finales','Fin de las estimaciones','help'));

close all
 

