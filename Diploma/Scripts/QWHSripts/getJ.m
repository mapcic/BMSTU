% getJ: function description dx, meff, Ec, Ez
function J = getJ(dx, meff, Ec, dU, EFermi)
	e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); 
	hbar = 1.054*1e-34; k_B = 1.38e-23;
	T = 300;
	kT = T*k_B;

	k = ((2*meff(1)*e*kT)/(4*pi^2*hbar^3));
	J = k*ones(1, length(dU));

	for j = 1:length(dU)
		Uj = Ec - linspace( 0, dU(j), length(Ec) );
		dTDEz = @(Ez) TDEz(dx, meff, Uj, Ez, EFermi);
		J(j) = J(j)*integral(dTDEz, 0, max(1)*eVtoJ, 'AbsTol', 1e-36);
	end
end