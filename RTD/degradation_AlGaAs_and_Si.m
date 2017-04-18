%% degradation_AlGaAs_and_Si: function description
function [Ec, meff] = degradation_AlGaAs_and_Si(x_Al, C_Si, checkTime, dx, T)
	e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1);
	nm = 1e-9; me = 9.1*1e-31;
	hbar = 1.054*1e-34; k_B = 1.38e-23;

	kT = T*k_B; % J

	Time = max(checkTime)*365*24; % to hours

	n_Atoms = 4.42*1e28; % number Atoms in GaAs ~ AlAs
	n_Al = n_Atoms/2; % number atoms of Al in AlAs

	dt = 1; % one hour
	dtdx2 = dt*60*60/dx^2; % s/m^2

	D = 0.2*exp(-3.5/(kT*JtoEv))*1e-4; % m^2/s

	C_Al = x_Al*n_Al;
	len_C_Al = length(C_Al);

	d1 = D*dtdx2*ones(1, len_C_Al-1);
	d2 = [ 1 - D*dtdx2, 1 - 2*D*dtdx2*ones(1, len_C_Al-2), 1 - D*dtdx2 ];
	d3 = D*dtdx2*ones(1, len_C_Al-1);

	Matrix_Al = diag(d1, -1) + diag(d2) + diag(d3, +1);

	checkTime = [1, 5, 10, 15, 20, 25]*365*24;

	C_Al = C_Al';
	C_Al = C_Si';
	for j = 0 : dt : Time
		C_Al = Matrix_Al*C_Al;
		C_Si = Matrix_Si*C_Si;
		ind = find(j == checkTime); 
		if ind
			[Ec(ind, :), meff(ind, :)] = getEcAlGaAs(C_Al');
		end
	end