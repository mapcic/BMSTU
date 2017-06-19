function val = TDEz(delta, meff, U, Ez, EFermi)
	hbar = 1.0551*1e-34;
	k_B = 1.38e-23;
	T = 300;
	kT = T*k_B;

	kLeft = abs( sqrt( 2*meff(1)*(Ez-U(1)) )/hbar );
	kRight = abs( sqrt( 2*meff(end)*(Ez-U(end)) )/hbar );

	[waveLeft, ~] = getWaveFunction(delta, meff, U, Ez);

	T = (kRight./kLeft).*(meff(1)/meff(end)).*(abs(waveLeft(:, end)).^2)';
	D = log( ( 1 + exp( (EFermi + U(1) - Ez)/kT ) ) ./ ( 1 + exp( (EFermi + U(end) - Ez)/kT ) ) );

	val = T.*D;
end