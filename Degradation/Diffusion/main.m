clear; clc;

e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); 
me = 9.11*1e-31; nm = 1e-9;

T = 920; % K

checkTime = 0 : 5; %years

dx = 1; %nm

a = 10; % monolayers
b = 30; 

sizeHS = a + b + a;

grid_x_Al = [
	zeros(1, a), ...
	ones(1, b), ...
	zeros(1, a)
];

nSi = 1e24;

[~, ~, grids_C_Al] = getDiffCloseAlGaAs( grid_x_Al, checkTime, dx*nm, T );
showResult(dx, sizeHS, grids_C_Al, 'Close Diffusion System', checkTime);

[~, ~, grids_C_Al] = getDiffOpenAlGaAs( grid_x_Al, checkTime, dx*nm, T );
showResult(dx, sizeHS, grids_C_Al, 'Open Diffusion System', checkTime);

[~, ~, grids_C_Al, grids_C_Si] = getDiffAlGaAs_Si(grid_x_Al, checkTime, dx*nm, 430, nSi);
showResultSi(dx, sizeHS, grids_C_Al, grids_C_Si, 'Not Constant $D$', checkTime);