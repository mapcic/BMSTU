clear; clc;
// Standart constants
nm = 1e-9;
eVtoJ = 1.6e-19;
JtoEv = eVtoJ^(-1);
me = 9.10938356*1e-31;
hbar = 1.0551*1e-34;
// Grid step
delta = 0.3; // nm ~ atom's radius
// Count layers
a = 33;
b = 6;
c = 9;
// Effective mass
m_AlAS = 0.15;
m_GaAs = 0.067;
// Hight of potential barrier
U = 1;
// Applyied voltage
dU = 0;
// Potential grid
gridU = [zeros(1, a), U*ones(1, b), zeros(1, c), U*ones(1, b), zeros(1, a)];
// Mass grid
gridMeff = [m_GaAs*ones(1, a), m_AlAS*ones(1, b), m_GaAs*ones(1, c), m_AlAS*ones(1, b), m_GaAs*ones(1, a)];
// Energy
dEz = 0.1;
EzRange = [0, U];
// Energy Fermi