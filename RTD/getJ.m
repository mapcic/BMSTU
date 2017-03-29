% getJ: function description delta, meff, U, Ez
function J = getJ(delta, meff, U, dU, EFermi)
	e = 1.6e-19;
	hbar = 1.054*1e-34;
	k_B = 1.38e-23;
	T = 300;
	kT = T*k_B;

	k = ((2*meff(1)*e*kT)/(4*pi^2*hbar^3));
	J = k*ones(1, length(dU));

	for j = 1:length(dU)
		Uj = U - linspace( 0, dU(j), length(U) );
		dTDEz = @(Ez) TDEz(delta, meff, Uj, Ez, EFermi);
		J(j) = J(j)*integral(dTDEz, 0, max(Uj), 'AbsTol', 1e-30);
	end
end