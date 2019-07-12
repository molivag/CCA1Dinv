function [fk] = ker1(r,e,n,nlam,rab,imi)
s = 1.7239458;
n1=n-1;
ss=log(10)/10;
%ss = 1;
abnor=rab; 
for  i=1:nlam
      exmod=i - imi + log(10*abnor)/ss;
      al=(1/s)*exp(-ss*exmod);
      e1=-2.*al*e(n1);
      tn=r(n);
      for j=n:-1:2 
      rhoj=r(j-1);
      ta1=tanh(al*e(j-1));
      tn=(tn + rhoj*ta1)/(1 + tn*ta1/rhoj);
      end
      fk(i)=tn;
      end
