clear; clc;

e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); 
me = 9.11*1e-31; nm = 1e-9;
hbar = 1.054*1e-34; k_B = 1.38e-23;
T = 300;
kT = T*k_B;

% Do it smooth
dis = 1;
% atoms' radius
dx = 0.56; %nm
dx = dx/dis;
 
L = 20;
R = 20;
HS = L + R;

% Fermi Energy
EFermi = 1.51*1e-20; % J

[Ec, Eg, me_eff, mp_eff] = getBandPropAlGaAs( [ones(1, L), zeros(1, R)] );

ni = 1e12;
Nd = 1e22;
Na = 1e20;

Ndi = [Nd*ones(1, L), ni*ones(1, R)];
Nai = [ni*ones(1, L), Na*ones(1, R)];

eps = 13.18 - 3.12*[ones(1, L), zeros(1, R)];

Nc = 2*(me*me_eff*kT/pi/hbar^2/2).^(3/2);
Nv = 2*(me*mp_eff*kT/pi/hbar^2/2).^(3/2);

ni = sqrt(Nc.*Nv).*exp(-Eg./(2*kT*JtoEv));

E_F_n = kT*log(Nd./ni(1:L))*JtoEv;
E_F_p = - kT*log(Na./ni(L+1:R+L))*JtoEv;

% E_F_nd = (3*pi^2*Nd*ones(1, L)).^(2/3)*hbar^(2)./(2*me*me_eff(1:L))*JtoEv;
% E_F_pd = (3*pi^2*Na*ones(1, R)).^(2/3)*hbar^(2)./(2*me*mp_eff(L+1:R+L))*JtoEv;

subplot(1, 4, 1);
plot(0:L-1, Eg(1:L)/2,...
	L:L+R-1, Eg(L+1:R+L)/2,...
	0:L-1, -Eg(1:L)/2,...
	L:L+R-1, -Eg(L+1:R+L)/2,...
	0:L-1, E_F_n,...
	L:L+R-1, E_F_p);

subplot(1, 4, 2);
plot(0:L-1, Eg(1:L)/2 - E_F_n,...
	L:L+R-1, Eg(L+1:R+L)/2 - E_F_p,...
	0:L-1, -Eg(1:L)/2 - E_F_n,...
	L:L+R-1, -Eg(L+1:R+L)/2 - E_F_p);

Ec = [Eg(1:L)/2 - E_F_n, Eg(L+1:R+L)/2 - E_F_p];
Ev = [-Eg(1:L)/2 - E_F_n, -Eg(L+1:R+L)/2 - E_F_p];

subplot(1, 4, 3);
plot(0:HS-1, Ec,...
	0:HS-1, Ev);

% [V, n] = getConcentrationElectrons(0.001,...
% 	Ec,...
% 	meff,...
% 	Ndi,...
% 	eps,...
% 	dx,...
% 	L...
% );

% Uj = Ec - V*eVtoJ;

% showResult(dx*nm, sizeHS, Ec1, Ec2, Ec3, Ec4, J1, J2, J3, J4, dU, Tr1, Tr2, Tr3, Tr4, a);
% showResult(dx*nm, sizeHS, Ec1, Ec2, Ec3, Ec4, J1, J2, J3, J4, dU, a);
