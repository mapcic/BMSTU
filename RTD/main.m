clear; clc;

e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); 
me = 9.11*1e-31; nm = 1e-9;

T = 800; % K

checkTime = [1, 5, 10, 15, 20, 25]; %years

% atoms' radius
dx = 0.56; %nm

% Count layers
a = 8; % monolayers
b = 5;
c = 6;

% Fermi Energy
EFermi = 1.51*1e-20; % J

% Applyied voltage
dU = 0:0.01:0.6;

% grid of Al conentration
grid_x_Al = [zeros(1, a), ones(1, b), zeros(1, c), ones(1, b), zeros(1, a)];

% Get profile Ec
[grids_Ec, grids_meff] = degradation_AlGaAs( grid_x_Al, checkTime, dx*nm, T );

% get J from V
for j = 1 : length(checkTime)
	J = getJ(dx*nm, grids_meff(j, :)*me, grids_Ec(j, :)*eVtoJ, dU*eVtoJ, EFermi);
	plot(dU, J);
	hold on;
end