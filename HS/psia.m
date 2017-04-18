function [psiL,psiR]=psia(Ui,mi,delta,Ez)
	hp=1.054E-34;
	n=length(Ui);
	psiL=zeros(length(Ez),n);
	psiR=zeros(length(Ez),n);
	for I=1:length(Ez)
	kl=sqrt(2*mi(1)*(Ez(I)-Ui(1)))/hp;
	kr=sqrt(2*mi(end)*(Ez(I)-Ui(end)))/hp;
	c1=ones(n-1,1);
	c2=((2*(delta^2)*mi(2:end-1)./(hp^2)).*(Ez(I)-Ui(2:end-1)))-(mi(2:end-1)./mi(3:end))-1;
	c2=[1i*kl*delta-1,c2,1i*kr*delta-1];
	c3=[1,mi(2:end-1)./mi(3:end)];
	H=diag(c1,-1)+diag(c2)+diag(c3,1);
	fl=[2*1i*kl*delta;zeros(n-1,1)];
	fr=[zeros(n-1,1);2*1i*kr*delta];
	psiL(I,:)=((H^(-1))*fl)';
	psiR(I,:)=((H^(-1))*fr)';
	end
end