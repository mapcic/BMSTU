clear; clc;

e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); 
me = 9.11*1e-31; nm = 1e-9;

T = 100 + 273.5; % K

checkTime = [0, 1, 5, 15]*360*24; %seconds

% atoms' radius
smth = 4;
dx = 0.56/smth; %nm

% Count layers
a = 10*smth; % monolayers
b = 6*smth;
c = 6*smth;

sizeHS = a + b + c + b + a;


% Fermi Energy
EFermi = 1.51*1e-20; % J

% Applyied voltage
dU = 0:0.25:0.5;

% grid of Al conentration
grid_x_Al = [zeros(1, a), ...
	0.44*ones(1, b), ...
	zeros(1, c), ...
	0.44*ones(1, b), ...
	zeros(1, a)
];

% Get profile Ec
[grids_Ec, grids_meff, grids_C_Al] = getDiffBOpenAlGaAs( grid_x_Al, checkTime, dx*nm, T );

% get J from V
for j = 1 : length(checkTime)
	J(j, :) = getJ(dx*nm, ...
		grids_meff(j, :)*me, ...
		grids_Ec(j, :)*eVtoJ, ...
		dU*eVtoJ, ...
		EFermi...
	);
end

showResult(grids_C_Al, grids_Ec, checkTime, J, sizeHS, dx, T, dU); % Img 4.3, Img 4.4