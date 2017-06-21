% img 2.5, 2.7, 2.9
clear; clc;

e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); 
me = 9.11*1e-31; nm = 1e-9;
hbar = 1.054*1e-34; k_B = 1.38e-23;

TmFrid = [0, 4, 10, 18, 28]; %Times grid. When check diffusion result?

T = 273+300; % Temperature in Kelvin
kT = k_B*T;

dx = 0.5; %nm

a = 200; % monolayers
b = 60;
c = 3;
d = 3;
e = 9;
f = 6;
g = 20;
k = 60;
m = 200; 

HS = a + b + c + d + e + f + g + k + m;

nAtoms = 4.2*1e28; % number Atoms in GaAs ~ AlAs
nAl = nAtoms/2; % number atoms of Al in AlAs

EgGaAs = 1.521 - 5.58*1e-4*T^2/(220+T);
EgAlAs = 2.239 - 6*1e-4*T^2/(T+408);

NcGaAs = 2*(me*0.067*kT/pi/hbar^2/2)^(3/2);
NvGaAs = 2*(me*0.49*kT/pi/hbar^2/2)^(3/2);

NcAlAs = 2*(me*0.15*kT/pi/hbar^2/2)^(3/2);
NvAlAs = 2*(me*0.81*kT/pi/hbar^2/2)^(3/2);

niGaAs = sqrt(NcGaAs*NvGaAs)*exp(-EgGaAs/(2*kT*JtoEv));
niAlAs = sqrt(NcAlAs*NvAlAs)*exp(-EgAlAs/(2*kT*JtoEv));

niGrid = [
	niGaAs*ones(1, a), ...
	niGaAs*ones(1, b), ...
	niGaAs*ones(1, c), ...
	niAlAs*ones(1, d), ...
	niGaAs*ones(1, e), ...
	niAlAs*ones(1, f), ...
	niGaAs*ones(1, g), ...
	niGaAs*ones(1, k), ...
	niGaAs*ones(1, m)
];

AlGrid = [
	zeros(1, a), ...
	zeros(1, b), ...
	zeros(1, c), ...
	ones(1, d), ...
	zeros(1, e), ...
	ones(1, f), ...
	zeros(1, g), ...
	zeros(1, k), ...
	zeros(1, m)
];
AlGrid = nAl*AlGrid;

nSip = 5e24;
nSi = 2e23;
SiGrid = [
	nSip*ones(1, a), ...
	nSi*ones(1, b), ...
	niGaAs*ones(1, c), ...
	niAlAs*ones(1, d), ...
	niGaAs*ones(1, e), ...
	niAlAs*ones(1, f), ...
	niGaAs*ones(1, g), ...
	nSi*ones(1, k), ...
	nSip*ones(1, m)
];

[diffAlGrid, diffSiGrid] = getDiffusion(AlGrid, SiGrid, niGrid, TmGrid, dx*nm, T);
% showResultSi(dx, sizeHS, grids_C_Al, grids_C_Si, 'Not Constant $D$', checkTime);