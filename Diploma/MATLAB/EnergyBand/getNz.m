% Fermi level is 0eV.
function [nz, pz] = getNz(Un, Up, meff_e, meff_p)
	e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); 
	hbar = 1.0551*1e-34; k_B = 1.38e-23;
	T=300;

	kT = k_B*T;

	Nc3D = 4*pi*(2*meff_e./(2*pi*hbar)^2)^(3/2);
	Nv3D = 4*pi*(2*meff_p./(2*pi*hbar)^2)^(3/2);

	for j = 1 : length(U)
		foo = @(Ez) sqrt(Ez - U(j))./(1 + exp((Ez)/(kT)));
		nz(j) = Nc3D(j)*integral(foo, Un(j), Un(j) + 2*eVtoJ, 'AbsTol', 1e-100);

		foo = @(Ez) sqrt(Ez - U(j))./(1 + exp((Ez)/(kT)));
		pz(j) = Nv3D(j)*integral(foo, Up(j), Up(j) + 2*eVtoJ, 'AbsTol', 1e-100);
	end
end