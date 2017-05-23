function Vnew = solvePoisonEq(dU, Vold, nold, eps, Ni, dx)
	e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); 	
	k_B = 1.38e-23; eps0 = 8.85e-12;

	T = 300;
	Vref = k_B*T*JtoEv;

	lenV = length(Vold);
	
	d1 = [ones(1, lenV-2), 0];
	
	d2 = ...
		-1 ...
		-(eps(1:end-2)./eps(3:end))...
		-e*dx^2/(eps0*Vref)./eps(3:end).*nold(2:end-1);
	d2 = [1, d2, 1];

	d3 = eps(1:end-2)./eps(3:end);
	d3 = [0, d3];
	
	dfree = ...
		e/eps0*dx^2 ...
		./eps(1:end-2)...
		.*(nold(2:end-1).*(1-Vold(2:end-1)/Vref)-Ni(2:end-1));
	dfree = [0, dfree, dU*JtoEv]';
	
	Matrix = diag(d1, -1) + diag(d2) + diag(d3, 1);
	Vnew = (Matrix\dfree)';
end
