%% getDODS_AlGaAs: function description
function [Ec, meff, Alx] = getDODS_AlGaAs(dx, checkTime, T, hsSize, hsStep, hsAl)
	e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1);
	nm = 1e-9; me = 9.1*1e-31;
	hbar = 1.054*1e-34; k_B = 1.38e-23;

	kT = T*k_B; % J

	D_Al = 0.2*exp(-3.5/(kT*JtoEv))*1e-4; % m^2/s

	Time = max(checkTime)*12*30*24; % to hours

	n_Atoms = 4.42*1e28; % number Atoms in GaAs ~ AlAs
	n_Al = n_Atoms/2; % number atoms of Al in AlAs

	C_Al = [];
	dxGrid = [];

	for ind = 1 : length(hsSize)
		tmp = ones( 1, hsSize(ind)/hsStep(ind));
		C_Al = [C_Al, hsAl(ind)*tmp];
		dxGrid = [dxGrid, dx*hsStep(ind)*tmp];
	end

	C_Al = C_Al*n_Al;
	% dxGrid = dxGrid(1:end-1);
	len = length(C_Al);

	dt = 60*60; % one hour
	% dxpm = [ dxGrid(1), (dxGrid(2:end-1) + dxGrid(1:end-2)) / 2, dxGrid(end)];
	dxpm = (dxGrid(2:end-1) + dxGrid(1:end-2)) / 2;
	
	% dxGrid
	% dxpm

	lmbd = D_Al*dt./dxpm; 
	d1 = [lmbd./dxGrid(1:end-2), 0];
	d2 = [ 1, 1 - lmbd.*(1./dxGrid(1:end-2) + 1./dxGrid(2:end-1)), 1 ];
	d3 = [0, lmbd./dxGrid(2:end-1)];

	Matrix_Al = diag(d1, -1) + diag(d2) + diag(d3, +1);

	if (find(0 == checkTime))
		[Ec(1, :), ~, meff(1, :), ~] = getBandPropAlGaAs(C_Al);
		Alx(1, :) = C_Al./n_Al;		
	end
	% Ec = 1;
	% meff = 1;
	% Alx = 1;

	C_Al = C_Al';
	for j = 1 : Time
		% clc; disp(j/Time*100);
		C_Al = Matrix_Al*C_Al;
		ind = find(j == checkTime*12*30*24); 
		if (ind & j ~= 0)
			[Ec(ind, :), ~, meff(ind, :), ~] = getBandPropAlGaAs(C_Al');
			Alx(ind, :) = C_Al'./n_Al;
		end
	end
end