%% getDiffOpenAlGaAsNd: function description
function [resAlGrid, resSiGrid] = getDiffAlGaAs_Si(AlGrid, SiGrid, niGrid, TmGrid, dx, T)
	e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1);
	nm = 1e-9; me = 9.1*1e-31;
	hbar = 1.054*1e-34; k_B = 1.38e-23;

	kT = T*k_B; % J

	time = max(TmGrid); % to
	time

	dt = 1; % step grid of time
	% dtdx2 = dt*60*60/dx^2; % s/m^2, T is hours
	% dtdx2 = dt*60/dx^2; % s/m^2, T is minute
	dtdx2 = (dt*1e-5)/dx^2; % s/m^2, T is sendond

	if (find(0 == TmGrid))
		resAlGrid(1, :) = AlGrid;		
		resSiGrid(1, :) = SiGrid;		
	end

	AlGrid = AlGrid'; % To multiply
	SiGrid = SiGrid';

	for j = 0 : dt : time
		D = 0.17*exp(-3.5/(kT*JtoEv))*(SiGrid'./niGrid).^3;

		d1 = [D(1:end-2)*dtdx2, 0];
		d2 = [1, 1 - 2*D(2:end-1)*dtdx2, 1];
		d3 = [0, D(2:end-1)*dtdx2];
		Matrix = diag(d1, -1) + diag(d2) + diag(d3, +1);

		AlGrid = Matrix*AlGrid;
		SiGrid = Matrix*SiGrid;

		ind = find(j == TmGrid); 
		if (ind & j ~= 0)
			resAlGrid(ind, :) = AlGrid;
			resSiGrid(ind, :) = SiGrid;
		end
	end
end