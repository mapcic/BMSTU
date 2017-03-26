clear; clc;

delta = 0.56; % nm ~ atom's radius
% Count layers
a = 8;
b = 4;
c = 6;
% Effective mass
m_AlAS = 0.15;
m_GaAs = 0.067;
% Hight of potential barrier
U = 1;
% Applyied voltage
dU = 0:0.001:0.6;
%  Potential grid
gridU = [zeros(1, a), U*ones(1, b), zeros(1, c), U*ones(1, b), zeros(1, a)];
% Mass grid
gridMeff = [m_GaAs*ones(1, a), m_AlAS*ones(1, b), m_GaAs*ones(1, c), m_AlAS*ones(1, b), m_GaAs*ones(1, a)];

J = getJ(delta, gridMeff, gridU, dU, 0.1);

plot(dU, J);