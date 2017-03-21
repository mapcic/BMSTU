%% getJ: function description delta, meff, U, Ez
function [J] = getJ(delta, meff, U, dU, EFermi)
	e = 1.6e-19;
	eVtoJ = e;
	JtoEv = e^(-1);
	me = 9.10938356*1e-31;
	hbar = 1.0551*1e-34;
	nm = 1e-9;
	T = 300;
	k_B = 1.38e-23;
	kT = T*k_B;

	EFermi = EFermi*eVtoJ;


	function val = TDEz(delta, meff, U, Ez)
		kLeft = abs( sqrt( 2*meff(1)*me*(Ez(j)-U(1)) )/hbar );
		kRight = abs( sqrt( 2*meff(end)*me*(Ez(j)-U(1)) )/hbar );

		[waveLeft, ~] = getWaveFunction(delta, meff, U, Ez);

		T = (kRight/kLeft)*(mi(1)/mi(end))*(abs(waveLeft(:,end)).^2)';
		D = log( 
				( 1 + exp( (EFermi + U(1) - Ez)/kT ) ) 
				./ 
				( 1 + exp( (EFermi + U(end) - Ez)/kT ) ) 
			);

		val = T.*D;
	end


	k = ((2*meff*me(1)*e*kT)/(4*pi^2*(hbar)^3));
	J = k*ones(1, length(dU));

	for j = 1:length(dU)
		dTDEz = @(Ez)TDEz(
				U - linspace( 0, dU(j),length(U) ),
				m, delta, Ez
			);
		J(j) = J(j)*integral(dTDEz, 0, 1*qe, 'AbsTol', 1e-30);
	end
end