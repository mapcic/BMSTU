% Img 4.5, 4.6
clear; clc;

e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); 
me = 9.11*1e-31; nm = 1e-9;

T = 100 + 273.5; % K
checkTime = [0, 1, 5, 15]*360*24; %years

% T = 800 + 273.5; % K
% checkTime = [0, 15, 30, 50]; %seconds

% T = 650 + 273.5; % K
% checkTime = [0, 560, 1200, 7200]; %seconds

smth = 4;
% atoms' radius
dx = 0.56/smth; %nm

% Count layers
a = 10*smth; % monolayers
b = 6*smth;
c = 6*smth;

sizeHS = a + b + c + b + a;

% grid of Al conentration
grid_x_Al = [zeros(1, a), ...
	0.44*ones(1, b), ...
	zeros(1, c), ...
	0.44*ones(1, b), ...
	zeros(1, a)
];

Nd = 5;
% Get profile Ec
[grids_Ec, grids_meff, grids_C_Al] = getDiffBOpenAlGaAsNd( grid_x_Al, checkTime, dx*nm, T, Nd );

showResult(grids_C_Al, grids_Ec, checkTime, sizeHS, dx, T, Nd);