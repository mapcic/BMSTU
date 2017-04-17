clear; clc;

e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); me = 9.11*1e-31;
nm = 1e-9; 
hbar = 1.054*1e-34; k_B = 1.38e-23;

T = 300;
kT = T*k_B;

me_eff = 0.062;
mp_eff = 0.5;

Nc = 2*(me*me_eff*kT/pi/hbar^2/2)^(3/2)*1e-6;
Nv = 2*(me*mp_eff*kT/pi/hbar^2/2)^(3/2)*1e-6;

Eg = 1.42;

ni = sqrt(Nc*Nv).*exp(-Eg./(2*kT*JtoEv));

Nd = 1e30;
Ef = 0.5*Eg + 1/2*kT*JtoEv*log(Nd/Nc);

n = Nc*exp(-Ef/kT/JtoEv)

grid_n_conduct = 0.5*1e6*( 2 + 0.5*(2*ni./(1e6)).^2 );