% calc integral by param

clear; clc;

L1 = 1e-5;
L2 = 8e-5;
t = 1000;

syms tau;

P = @(t) exp(-2*L1*t) + 2*exp(-L2*t).*double(int(L1*exp(tau*(L2 - 2*L1)), 0, t));
% 
P1 = double(P(t));
% 
T = integral(P, 0, inf);
