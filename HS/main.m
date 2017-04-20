clear; clc;

e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); 
me = 9.11*1e-31; nm = 1e-9;

T = 300; % K

checkTime = 0:0.5:3; %years

% atoms' radius
dx = 0.56; %nm

% Count layers
rez = 5;
a = 8; % monolayers
b = 5;
c = 6;

% Fermi Energy
EFermi = 1.51*1e-20; % J

% Applyied voltage
dU = 0:0.01:0.6;

Ec = [zeros(1, rez), zeros(1, a), ones(1, b), zeros(1, c), ones(1, b), zeros(1, a), zeros(1, rez)];
meff = [0.067*ones(1, rez), 0.067*ones(1, a), 0.15*ones(1, b), 0.067*ones(1, c), 0.15*ones(1, b), 0.067*ones(1, a), 0.067*ones(1, rez)];

Nd = 2.5*1e24;
ni = 1e12;
Ni = [Nd*ones(1, rez), ni*ones(1, a), ni*ones(1, b), ni*ones(1, c), ni*ones(1, b), ni*ones(1, a), Nd*ones(1, rez)];
accur = 0.000001;

eps = 13.18 - 3.12*Ec;

[V, n] = getConcentrationElectrons(accur, Ec*eVtoJ, meff*me, Ni, eps, dx*nm, [0], rez+1, length(Ec)-rez);
% [V, n] = UN(accur, Ec*eVtoJ, meff*me, Ni, eps, dx*nm, [0], rez+1, length(Ec)-rez);