clear; clc;

	e = 1.6e-19; eV2J = e; J2eV = e^(-1);
	hbar = 1.0551*1e-34; k_B = 1.38e-23;
	T = 300;

	meff = 0.15*me;
	Nd = 4*pi*(2*meff/(2*pi*hbar)^2)^(3/2);

	n = 
