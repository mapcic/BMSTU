%% getDiffOpenAlGaAsNd: function description
function [Ec, meff, Alx] = getDiffOpenAlGaAs(x_Al, checkTime, dx, T, Nd)
	e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1);
	nm = 1e-9; me = 9.1*1e-31;
	hbar = 1.054*1e-34; k_B = 1.38e-23;

	kT = T*k_B; % J

	Time = max(checkTime)*12*30*24; % to hours

	n_Atoms = 4.42*1e28; % number Atoms in GaAs ~ AlAs
	n_Al = n_Atoms/2; % number atoms of Al in AlAs

	dt = 1; % one hour
	dtdx2 = dt*60*60/dx^2; % s/m^2

	Eg_GaAs = 1.519 - 5.405*1e-4*T^2/(T+204);
	Nc = 2*(me*0.067*kT/pi/hbar^2/2)^(3/2);
	Nv = 2*(me*0.51*kT/pi/hbar^2/2)^(3/2);
	ni = sqrt(Nc*Nv)*exp(-Eg_GaAs/(2*kT*JtoEv));

	D_Al = 0.2*exp(-3.5/(kT*JtoEv))*1e-4*(Nd/ni)^3; % m^2/s

	C_Al = x_Al*n_Al;
	len = length(x_Al);

	d1 = [D_Al*dtdx2*ones(1, len-2), 0];
	d2 = [ 1, 1 - 2*D_Al*dtdx2*ones(1, len-2), 1 ];
	d3 = [0, D_Al*dtdx2*ones(1, len-2)];
	
	Matrix_Al = diag(d1, -1) + diag(d2) + diag(d3, +1);

	if (find(0 == checkTime))
		[Ec(1, :), ~, meff(1, :), ~] = getBandPropAlGaAs(C_Al);
		Alx(1, :) = C_Al./n_Al;		
	end

	C_Al = C_Al';
	for j = 0 : dt : Time
		% clc; disp(j/Time*100);
		C_Al = Matrix_Al*C_Al;
		ind = find(j == checkTime*12*30*24); 
		if (ind & j ~= 0)
			[Ec(ind, :), ~, meff(ind, :), ~] = getBandPropAlGaAs(C_Al');
			Alx(ind, :) = C_Al'./n_Al;
		end
	end
end