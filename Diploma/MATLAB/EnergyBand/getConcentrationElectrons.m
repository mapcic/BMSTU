function [Vnew, nold] = getConcentrationElectrons(accur, Ec, meff, Ni, eps, dx)
	e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); 	

	lenEc = length(Ec);
	
	Vnew = zeros(1, lenEc);
	Vold = Vnew + 10;
	
	while ( max(abs(Vnew - Vold)) > 1e-100 )
		Vold = Vnew;

		Ui = Ec - Vold*eVtoJ;

		nold = getNz(Ui, meff);

    	Vnew = solvePoisonEq(Vold, nold, eps, Ni, dx);
	end
	Ui = Ec - Vnew*eVtoJ;
	nold = getNz(Ui, meff);
end