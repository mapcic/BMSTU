function nz = getNz(Ui, meff, dx, boundL, boundR)
	e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); 
	hbar = 1.0551*1e-34; k_B = 1.38e-23;
	T=300;

	EFermi=1.51e-20;

	Nc3D = 4*pi*(2*meff(1)/(2*pi*hbar)^2)^(3/2);
	Nc3DActive = sqrt(2)*meff(end)^(3/2)*k_B*T/((2*pi)^2*hbar^3);
	
	function res = NEz(Ui, meff, dx, Ez, U1, Un)
		[waveL, waveR] = getWaveFunction(dx, meff, Ui, Ez);

		waveL = abs(waveL).^2;
		waveR = abs(waveR).^2;
		
		Ez = repmat(Ez', 1, length(Ui));
		
		waveL = Nc3DActive*(waveL)./sqrt(Ez - Ui(1)).*log(1 + exp((EFermi + U1 - Ez)/(k_B*T)));
		waveL( Ez <= Ui(1) ) = 0;
		
		waveR = Nc3DActive*waveR./sqrt(Ez - Ui(end)).*log(1 + exp((EFermi + Un - Ez)/(k_B*T)));
		waveR( Ez <= Ui(end) ) = 0;
		
		res = waveL + waveR;
	end

	foo = @(Ez) NEz( Ui(boundL: boundR), meff(boundL: boundR), dx, Ez, Ui(1), Ui(end) );
	nzA = integral(foo, Ui(end), 2*e, 'AbsTol', 1E-100, 'ArrayValued', true);

	U1 = Ui(1: boundL-1);
	U2 = Ui(boundR + 1: end);

	nzL = zeros(1, length(U1));
	nzR = zeros(1, length(U2));

	for j = 1 : length(U1)
		foo = @(Ez) sqrt(Ez - U1(j))./(1 + exp((Ez - (EFermi + Ui(1)))/(k_B*T)));
		nzL(j) = Nc3D*integral(foo, U1(j), 2*e, 'AbsTol', 1e-200);
	end

	for j = 1 : length(U2)
		foo = @(Ez) sqrt(Ez - U2(j))./(1 + exp((Ez - (EFermi + Ui(end)))/(k_B*T)));
		nzR(j) = Nc3D*integral(foo, U2(j), 2*e, 'AbsTol', 1e-200);
	end

	nz = [nzL, nzA, nzR];
end
