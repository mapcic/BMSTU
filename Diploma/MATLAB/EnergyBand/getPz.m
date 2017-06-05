% Fermi level is 0eV.
function pz = getPz(Up, meff_p)
	e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); 
	hbar = 1.0551*1e-34; k_B = 1.38e-23;
	T=300;

	kT = k_B*T;

	Nv3D = 4*pi*(2*meff_p./(2*pi*hbar)^2)^(3/2);

	for j = 1 : length(Up)
		foo = @(Ez) sqrt( - Ez + Up(j) )./(1 + exp((-Ez)/(kT)));
		pz(j) = Nv3D(j)*integral(foo, Up(j) - 2*eVtoJ, Up(j), 'AbsTol', 1e-100);
	end
end