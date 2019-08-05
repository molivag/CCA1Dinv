
%Regularizado
 I=diag(ones(1,n));
 J_reg=J;
 Xmc_reg=rho';
 alfa=0.011;
 for i=1:3
     disp(['Iteraci???n: ', num2str(i),' XmcR']),disp(' ')
    pause(0.5)
    Resp_Xmc_reg=rhoa(ab,Xmc_reg,e,capas,m)';
    Xmc_reg= Xmc_reg + inv(J_reg'*J_reg + (alfa*I)) * J_reg' * (sev7 - Resp_Xmc_reg)
    Y_Xmc_reg = rhoa(ab,Xmc_reg,e,capas,m)'
    pause(0.5)
    J_reg= Jacobiano( Xmc_reg, h, Y_Xmc_reg, ab, e, capas, m );
    pause(0.5)
    loglog(ab,Y_Xmc_reg,'.-'),grid on, hold on, pause(2)
    stairs(z,Xmc_reg,'--')
    %rho = Xmc_reg';
     
 end
legend('Dat_{obs}','Y(x)_{0}','\rho_{x_{0}}','Y(x)_{reg1}','\rho_{x_{reg1}}',...
       'Y(x)_{reg2}','\rho_{x_{reg2}}','Y(x)_{reg3}','\rho_{x_{3}}')

%graficas finales Xmc_regularizado 
 figure('name','Ajuste Final Xmcp')
 errorbar(ab,sev7,des_sev7,'xk'),set(gca,'xscale','log'),set(gca,'yscale','log')
 hold on, grid on, xlabel('AB/2'),ylabel('\rho_{a}')
 loglog(ab,Resp,ab,Y_Xmc_reg)
 legend('Dat_{obs}','Y(x)_{0}','Y(x)_{pond}')
 figure('name','Modelo Final Xmc_regu')
 stairs(z,Xmc_reg,'--')
%Residuales
 figure('name','Residuales Xmc_reg')
 rang=[1:m];
 resi_Xmc_reg = sev7 - Y_Xmc_reg;, O=zeros(1,m);
 plot(rang, resi_Xmc_reg, rang, O,'k'),xlim([0 18])

uiwait(msgbox('Inicio de la Descomposicion en Valores Singulares de J, presione OK','Fin de Xmc regularizado','help'));

close all
