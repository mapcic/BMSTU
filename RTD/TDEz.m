function val = TDEz(delta, meff, U, Ez, EFermi)
	e = 1.6e-19;
	eVtoJ = e;
	JtoEv = e^(-1);
	me = 9.10938356*1e-31;
	hbar = 1.0551*1e-34;
	T = 300;
	k_B = 1.38e-23;
	kT = T*k_B;

	kLeft = abs( sqrt( 2*meff(1)*me*(Ez-U(1))*eVtoJ )/hbar );
	kRight = abs( sqrt( 2*meff(end)*me*(Ez-U(end))*eVtoJ )/hbar );

	[waveLeft, ~] = getWaveFunction(delta, meff, U, Ez);

	T = (kRight./kLeft).*(meff(1)/meff(end)).*(abs(waveLeft(:, end)).^2)';
	D = log( ( 1 + exp( eVtoJ*(EFermi + U(1) - Ez)/kT ) ) ./ ( 1 + exp( eVtoJ*(EFermi + U(end) - 1)/kT ) ) );

	val = T.*D;
end

% 	function res=TDEz(Ui,mi,delta,Ez)
% 		kl=sqrt(2*mi(1)*(Ez-Ui(1)))/hp;
% 		kr=sqrt(2*mi(end)*(Ez-Ui(end)))/hp;
		
% 		[psiL,~]=psia(Ui,mi,delta,Ez);
		
% 		TE=(abs(kr)./abs(kl)).*(mi(1)/mi(end)).*(abs(psiL(:,end)).^2)';
% 		DE=log((1+exp((Ef+Ui(1)-Ez)/(kb*T)))./(1+exp((Ef+Ui(end)-1)/(kb*T))));
% 		res=TE.*DE;
% 	end