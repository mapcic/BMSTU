clear; clc;

e = 1.6e-19;
eVtoJ = e;
JtoEv = e^(-1); 
me = 9.10938356*1e-31;
nm = 1e-9;

% atoms' radius
delta = 0.56;

% Count layers
a = 8;
b = 4;
c = 6;

% Effective mass
m_AlAS = 0.15;
m_GaAs = 0.067;

% Hight of potential barrier
U = 1;

% Fermi Energy
EFermi = 1.51*1e-20;

% Applyied voltage
dU = 0:0.01:0.5;

%  Potential grid
gridU = [zeros(1, a), U*ones(1, b), zeros(1, c), U*ones(1, b), zeros(1, a)];

% Mass grid
gridMeff = [m_GaAs*ones(1, a), m_AlAS*ones(1, b), m_GaAs*ones(1, c), m_AlAS*ones(1, b), m_GaAs*ones(1, a)];


% diffusion

% get J from V
% J = getJ(delta*nm, gridMeff*me, gridU*eVtoJ, dU*eVtoJ, EFermi);
% plot(dU, J);