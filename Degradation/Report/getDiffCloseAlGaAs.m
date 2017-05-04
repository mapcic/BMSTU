%% getDiffCloseAlGaAs: function description
function [Ec, meff] = getDiffCloseAlGaAs(dx, x_Al, T, checkTime)
	e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1);
	nm = 1e-9; me = 9.1*1e-31;
	hbar = 1.054*1e-34; k_B = 1.38e-23;

	kT = T*k_B; % J

	Time = max(checkTime)*12*30*24; % to hours

	n_Atoms = 4.42*1e28; % number Atoms in GaAs ~ AlAs
	n_Al = n_Atoms/2; % number atoms of Al in AlAs

	dt = 1; % one hour
	dtdx2 = dt*60*60/dx^2; % s/m^2

	D_Al = 0.2*exp(-3.5/(kT*JtoEv))*1e-4; % m^2/s

	x_Al = [0, zeros(1, a), 0.45*ones(1, b), zeros(1, c), 0.45*ones(1, b), zeros(1, a), 0];
	C_Al = x_Al*n_Al;
	len_C_Al = length(C_Al);

	d1 = D_Al*dtdx2*ones(1, len_C_Al-1);
	d2 = [ 1 - D_Al*dtdx2, 1 - 2*D_Al*dtdx2*ones(1, len_C_Al-2), 1 - D_Al*dtdx2 ];
	d3 = D_Al*dtdx2*ones(1, len_C_Al-1);
	Matrix_Al = diag(d1, -1) + diag(d2) + diag(d3, +1);

	checkTime = (1:25)*12*30*24;

	C_Al = C_Al';
	for j = 0 : dt : Time
		% clc; disp(j/Time*100);
		C_Al = Matrix_Al*C_Al;
		ind = find(j == checkTime); 
		if ind
			plot(1:len_C_Al, C_Al'/n_Al);
			hold on;
		end
	end
end