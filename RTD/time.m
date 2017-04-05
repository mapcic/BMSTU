clear; clc;

e = 1.6e-19;
eVtoJ = e;
JtoEv = e^(-1); 
nm = 1e-9;
me = 9.10938356*1e-31;
hbar = 1.054*1e-34;
k_B = 1.38e-23;
T = 600;
kT = T*k_B;

Eg = 1.519 - 5.405*1e-4*T^2/(T+204);

me_eff = 0.067;
mp_eff = 0.51;

Nc = 2*(me*me_eff*kT/pi/hbar^2/2)^(3/2);
Nv = 2*(me*mp_eff*kT/pi/hbar^2/2)^(3/2);

ni = sqrt(Nc*Nv)*exp(-Eg/(2*kT*JtoEv));

delta = 0.56*nm;

nSi = 5*1e18*1e6;

% D = 0.2*exp(-3.5/(kT*JtoEv))*(nSi/ni)^3;
D = 0.2*exp(-3.5/(kT*JtoEv));


dt = (delta)^2/(2*D*1e-4)/60/60/24;