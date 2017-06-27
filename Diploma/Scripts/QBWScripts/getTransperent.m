%% getTransperent: function description
function T = getTransperent(delta, meff, U, numPoint)
	e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); 
	hbar = 1.0551*1e-34;

	Ez = linspace(min(U), max(U), numPoint);
	% Ez = linspace(0.01*eVtoJ, max(U), numPoint);

	kLeft = abs( sqrt( 2*meff(1)*(Ez-U(1)) )/hbar );
	kRight = abs( sqrt( 2*meff(end)*(Ez-U(end)) )/hbar );

	[waveLeft, ~] = getWaveFunction(delta, meff, U, Ez);

	T = (kRight./kLeft).*(meff(1)/meff(end)).*(abs(waveLeft(:, end)).^2)';
end