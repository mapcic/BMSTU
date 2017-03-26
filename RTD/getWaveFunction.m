%% getWaveFunction: function description
function [waveLeft, waveRigth] = getWaveFunction(delta, meff, U, Ez)
	eVtoJ = 1.6e-19;
	JtoEv = eVtoJ^(-1);
	me = 9.10938356*1e-31;
	hbar = 1.0551*1e-34;
	nm = 1e-9;

	EzLen = length(Ez);
	ULen = length(U);

	Ez = Ez*eVtoJ;
	U = U*eVtoJ;
	delta = delta*nm;

	waveLeft = zeros(EzLen, ULen);
	waveRigth = zeros(EzLen, ULen);

	for j = 1 : EzLen
		kLeft = sqrt( 2*meff(1)*me*(Ez(j) - U(1)) )/hbar;
		kRight = sqrt( 2*meff(end)*me*(Ez(j) - U(end)) )/hbar;

		% d1 = meff(2:ULen)./meff(1:ULen-1);
		% d2 = 2*delta^2*me*meff(2:end-1).*(Ez(j) - U(2:end-1))./hbar^2 - 2 ;
		% d2=[1i*kLeft*delta - 1, d2, 1i*kRight*delta - 1];
		% d3 = meff(1:ULen-1)./meff(2:ULen);

		d1=ones(ULen-1,1);
		d2=((2*(delta^2)*me*meff(2:end-1)./(hbar^2)).*(Ez(j)-U(2:end-1)))-(meff(2:end-1)./meff(3:end))-1;
		d2=[1i*kLeft*delta-1,d2,1i*kRight*delta-1];
		d3=[1,meff(2:end-1)./meff(3:end)];

		H = diag(d1, -1) + diag(d2) + diag(d3, +1);

		fLeft = [2*1i*kLeft*delta; zeros(ULen-1, 1)];
		fRight = [zeros(ULen-1, 1); 2*1i*kRight*delta];

		waveLeft(j, :) = (inv(H)*fLeft)';
		waveRigth(j, :) = (inv(H)*fRight)';
	end
	
end


% function [psiL,psiR]=psia(Ui,mi,delta,Ez)
% 	hp=1.054E-34;
% 	n=length(Ui);
	
% 	psiL=zeros(length(Ez),n);
% 	psiR=zeros(length(Ez),n);
	
% 	for I=1:length(Ez)
% 		kl=sqrt(2*mi(1)*(Ez(I)-Ui(1)))/hp;
% 		kr=sqrt(2*mi(end)*(Ez(I)-Ui(end)))/hp;
		
% 		c1=ones(n-1,1);
% 		c2=((2*(delta^2)*mi(2:end-1)./(hp^2)).*(Ez(I)-Ui(2:end-1)))-(mi(2:end-1)./mi(3:end))-1;
% 		c2=[1i*kl*delta-1,c2,1i*kr*delta-1];
% 		c3=[1,mi(2:end-1)./mi(3:end)];
		
% 		H=diag(c1,-1)+diag(c2)+diag(c3,1);
		
% 		fl=[2*1i*kl*delta;zeros(n-1,1)];
% 		fr=[zeros(n-1,1);2*1i*kr*delta];
		
% 		psiL(I,:)=((H^(-1))*fl)';
% 		psiR(I,:)=((H^(-1))*fr)';
% 	end
% end