function [Vnew,nold] = getConcentrationElectrons(accur, Ec, meff, Ni, eps, dx, dU, boundL, boundR)
	e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); 	

	lenEc = length(Ec);
	
	Vnew = [zeros(1, boundL-1), linspace(0, dU, boundR - boundL + 1), dU*ones(1, lenEc - boundR)];
	Vold = Vnew + 10*eVtoJ;
	Vnew
	while ( max(abs(Vnew - Vold))*JtoEv > accur )
		Vold = Vnew;
		Ui = Ec - Vold;

		nold = getNz(Ui, meff, dx, boundL, boundR);
        
        Vnew = solvePoisonEq(dU, Vold, nold, eps, Ni, dx);
	end

	Ui = Ec - Vnew;
	nold = getNz(Ui, meff, dx, boundL, boundR);
end