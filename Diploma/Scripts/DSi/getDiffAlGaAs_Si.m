%% getDiffOpenAlGaAsNd: function description
function [Ec, meff, Alx, Six] = getDiffOpenAlGaAs(x_Al, checkTime, dx, T, n_Si)
	e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1);
	nm = 1e-9; me = 9.1*1e-31;
	hbar = 1.054*1e-34; k_B = 1.38e-23;

	kT = T*k_B; % J

	Time = max(checkTime); % to hours

	n_Atoms = 4.42*1e28; % number Atoms in GaAs ~ AlAs
	n_Al = n_Atoms/2; % number atoms of Al in AlAs

	dt = 1; % one hour
	dtdx2 = dt*60*60/dx^2; % s/m^2

	Eg_GaAs = 1.519 - 5.405*1e-4*T^2/(T+204);
	Nc = 2*(me*0.067*kT/pi/hbar^2/2)^(3/2);
	Nv = 2*(me*0.51*kT/pi/hbar^2/2)^(3/2);
	ni = sqrt(Nc*Nv)*exp(-Eg_GaAs/(2*kT*JtoEv));

	C_Al = [0, x_Al, 0]*n_Al;
	len = length(C_Al);

	C_Si = [n_Si, ni*ones(size(x_Al)), n_Si];

	if (find(0 == checkTime))
		[Ec(1, :), ~, meff(1, :), ~] = getBandPropAlGaAs(C_Al(2:end-1));
		Alx(1, :) = C_Al(2:end-1)./n_Al;		
		Six(1, :) = C_Si(2:end-1)./n_Si;		
	end

	C_Al = C_Al';
	C_Si = C_Si';
	for j = 0 : dt : Time
		D = 0.2*exp(-3.5/(kT*JtoEv))*(C_Si'./ni).^3*1e-4;

		D_plus = (D(1:len-1) + D(2:len))./2;
		D_minus = (D(2:len) + D(1:len-1))./2;

		d1 = [D_minus(1:end-1)*dtdx2, 0];
		d2 = [1, 1 - (D_plus(2:end) + D_minus(1:end-1))*dtdx2, 1];
		d3 = [0, D_plus(2:end)*dtdx2];
		Matrix = diag(d1, -1) + diag(d2) + diag(d3, +1);

		C_Al = Matrix*C_Al;
		C_Si = Matrix*C_Si;

		ind = find(j == checkTime); 
		if (ind & j ~= 0)
			[Ec(ind, :), ~, meff(ind, :), ~] = getBandPropAlGaAs(C_Al(2:end-1)');
			Alx(ind, :) = C_Al(2:end-1)'./n_Al;
			Six(ind, :) = C_Si(2:end-1)';
		end
	end
end