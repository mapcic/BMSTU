% img 2.5, 2.7, 2.9
clear; clc;

e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); 
me = 9.11*1e-31; nm = 1e-9;
hbar = 1.054*1e-34; k_B = 1.38e-23;

% TmGrid = [0, 4, 10, 18, 28]; %Times grid. When check diffusion result?
TmGrid = [0, 1]*60*60*1e5; %Times grid. When check diffusion result?

T = 600; % Temperature in Kelvin
kT = k_B*T;

dx = 0.5; %nm

a = 200; % monolayers
b = 60;

HS = a + b;

nAtoms = 4.2*1e28; % number Atoms in GaAs ~ AlAs
nAl = nAtoms/2; % number atoms of Al in AlAs

EgGaAs = 1.521 - 5.58*1e-4*T^2/(220+T);

NcGaAs = 2*(me*0.067*kT/pi/hbar^2/2)^(3/2);
NvGaAs = 2*(me*0.49*kT/pi/hbar^2/2)^(3/2);

niGaAs = sqrt(NcGaAs*NvGaAs)*exp(-EgGaAs/(2*kT*JtoEv));

niGrid = [
	niGaAs*ones(1, a), ...
	niGaAs*ones(1, b)
];

AlGrid = [
	zeros(1, a), ...
	zeros(1, b)
];
AlGrid = nAl*AlGrid;

nSip = 5e24;
nSi = 2e23;
SiGrid = [
	nSip*ones(1, a), ...
	nSi*ones(1, b)
];

[diffAlGrid, diffSiGrid] = getDiffusionSimple(AlGrid, SiGrid, niGrid, TmGrid, dx*nm, T);
% showResultSi(dx, sizeHS, grids_C_Al, grids_C_Si, 'Not Constant $D$', checkTime);