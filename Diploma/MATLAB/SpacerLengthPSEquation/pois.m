function Vnew=pois(V,Vold,nold,eps,Ni,delta)
T=300;
qe=1.602E-19;
eps0=8.85E-12;
kb=1.38E-23;
Vref=kb*T/qe;
n=length(Vold);
ai=[ones(1,n-2),0];
bi=-1-(eps(3:end)./eps(2:end-1))-(((delta^2)./eps(2:end-1).*(qe*nold(2:end-1)/(eps0*Vref))));
 bi=[1,bi,1];
 ci=eps(3:end)./eps(2:end-1);
 ci=[0,ci];
 di=((delta^2)./eps(2:end-1))*(qe/eps0).*(nold(2:end-1).*(1-Vold(2:end-1)/Vref)-Ni(2:end-1));
 di=[0,di,V]';
 C=diag(ai,-1)+diag(bi)+diag(ci,1);
 Vnew=C\di;
 Vnew=Vnew';
 end