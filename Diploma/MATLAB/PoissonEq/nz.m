function res=nz(Ui,mi,delta,i1,i2)
qe=1.602E-19;
m0=9.109E-31;
mw=0.067*m0;
Ef=1.51e-20;
kb=1.38E-23;
T=300;
hp=1.054E-34;
Cr=4*pi*((2*mw/((2*pi*hp)^2))^(3/2));
Ca=sqrt(2)*(mw^(3/2))*kb*T/(((2*pi)^2)*hp^3);
function res=NEz(Ui,mi,delta,Ez,U1,Un)
[pvL,pvR]=psia(Ui,mi,delta,Ez);
pvL=abs(pvL).^2;
pvR=abs(pvR).^2;
Ez=repmat(Ez',1,length(Ui));
pvL=Ca*((pvL)./sqrt(Ez-Ui(1))).*log(1+exp((Ef+U1-Ez)/(kb*T)));
pvL(Ez<=Ui(1))=0;
pvR=Ca*((pvR)./sqrt(Ez-Ui(end))).*log(1+exp((Ef+Un-Ez)/(kb*T)));
pvR(Ez<=Ui(end))=0;
res=pvL+pvR;
end
foo=@(Ez)NEz(Ui(i1:i2),mi(i1:i2),delta,Ez,Ui(1),Ui(end));
naz=integral(foo,Ui(end),2*qe,'AbsTol',1E-25,'ArrayValued',true);

U1=Ui(1:i1-1);
U2=Ui(i2+1:end);
nrz1=zeros(1,length(U1));
nrz2=zeros(1,length(U2));
for I=1:length(U1)
foo=@(Ez)sqrt(Ez-U1(I))./(1+exp((Ez-(Ef+Ui(1)))/(kb*T)));
nrz1(I)=Cr*integral(foo,U1(I),2*qe,'AbsTol',1E-25);
end
for I=1:length(U2)
foo=@(Ez)sqrt(Ez-U2(I))./(1+exp((Ez-(Ef+Ui(end)))/(kb*T)));
nrz2(I)=Cr*integral(foo,U2(I),2*qe,'AbsTol',1E-25);
end
res=[nrz1,naz,nrz2];
end