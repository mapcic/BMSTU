% getJ: function description delta, meff, U, Ez
function J = getJ(delta, meff, U, dU, EFermi)
	e = 1.6e-19;
	eVtoJ = e;
	JtoEv = e^(-1);
	me = 9.10938356*1e-31;
	hbar = 1.0551*1e-34;
	T = 300;
	k_B = 1.38e-23;
	kT = T*k_B;

	k = ((2*meff(1)*me*e*kT)/(4*pi^2*hbar^3));
	J = k*ones(1, length(dU));

	for j = 1:length(dU)
		Uj = U - linspace( 0, dU(j), length(U) );
		dTDEz = @(Ez) TDEz(delta, meff, Uj, Ez, EFermi);
		J(j) = e*J(j)*integral(dTDEz, 0, max(Uj), 'AbsTol', 1e-30);
	end
end


% function J=JVfd(V,Ec,mi,delta)
% 	Ef=1.51e-20;
% 	T=300;
% 	kb=1.38E-23;
% 	qe=1.602E-19;
% 	hp=1.054E-34;
	
% 	function res=TDEz(Ui,mi,delta,Ez)
% 		kl=sqrt(2*mi(1)*(Ez-Ui(1)))/hp;
% 		kr=sqrt(2*mi(end)*(Ez-Ui(end)))/hp;
		
% 		[psiL,~]=psia(Ui,mi,delta,Ez);
		
% 		TE=(abs(kr)./abs(kl)).*(mi(1)/mi(end)).*(abs(psiL(:,end)).^2)';
% 		DE=log((1+exp((Ef+Ui(1)-Ez)/(kb*T)))./(1+exp((Ef+Ui(end)-1)/(kb*T))));
% 		res=TE.*DE;
% 	end
	
% 	J=zeros(1,length(V));
	
% 	for I=1:length(V)
% 		Ui=Ec-qe*linspace(0,V(I),length(Ec));
		
% 		foo=@(Ez)TDEz(Ui,mi,delta,Ez);
		
% 		J(I)=((2*mi(1)*qe*kb*T)/(4*(pi^2)*(hp)^3))*integral(foo,0,1*qe,'AbsTol',1e-30);
% 	end
% end