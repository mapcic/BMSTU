clear; clc;

e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); 
me = 9.11*1e-31; nm = 1e-9;

T = 100 + 273.5; % K
checkTime = [0, 1, 5, 10, 11, 12, 13]*24*360; %years

smth = 4;
% atoms' radius
dx = 0.56/smth; %nm

% Count layers
r = 20*smth; % monolayers
a = 10*smth; % monolayers
b = 6*smth;
c = 6*smth;

sizeHS = r + a + b + c + b + a + r;

% grid of Al conentration
grid_x_Al = [zeros(1, a), ...
	0.44*ones(1, b), ...
	zeros(1, c), ...
	0.44*ones(1, b), ...
	zeros(1, a)
];

Ndp = 2.5*1e18*1e6;
Nd = 2*1e17*1e6;
% Get profile Ec
[~, ~, grids_C_Al, Six] = getDiffBAlGaAs_Si( grid_x_Al, checkTime, dx*nm, T, Ndp, Nd, r);

showResult(grids_C_Al, Six/Nd, grids_Ec, checkTime, sizeHS, dx, T, Nd);