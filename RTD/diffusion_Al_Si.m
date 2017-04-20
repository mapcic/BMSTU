clear; clc;

e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1);
nm = 1e-9; me = 9.1*1e-31;
hbar = 1.054*1e-34; k_B = 1.38e-23;

T = 350;
kT = T*k_B; % J

% Count layers
a = 8; % monolayers
b = 5;
c = 6;

Time = 15*12*30*24; % to hours

n_Atoms = 4.42*1e28; % number Atoms in GaAs ~ AlAs
n_Al = n_Atoms/2; % number atoms of Al in AlAs

Eg = 1.519 - 5.405*1e-4*T^2/(T + 204);

me_eff = 0.067;
mp_eff = 0.51;

Nc = 2*(me*me_eff*kT/pi/hbar^2/2)^(3/2);
Nv = 2*(me*mp_eff*kT/pi/hbar^2/2)^(3/2);

ni = sqrt(Nc*Nv)*exp(-Eg/(2*kT*JtoEv));

dx = 0.56*nm; %nm
dt = 1; % one hour
dtdx2 = dt*60*60/dx^2; % s/m^2

n_Si = 2*1e18*1e6;
C_Si = [n_Si, ni*ones(1, a), ni*ones(1, b), ni*ones(1, c), ni*ones(1, b), ni*ones(1, a), n_Si];

x_Al = [0, zeros(1, a), 0.45*ones(1, b), zeros(1, c), 0.45*ones(1, b), zeros(1, a), 0];
C_Al = x_Al*n_Al;
len_C_Al = length(C_Al);

checkTime = (0:3:15)*12*30*24;

C_Al = C_Al';
C_Si = C_Si';
for j = 0 : dt : Time
	D = 0.2*exp(-3.5/(kT*JtoEv))*(C_Si'./ni).^3*1e-4;

	D_plus = (D(1:len_C_Al-1) + D(2:len_C_Al))./2;
	D_minus = (D(2:len_C_Al) + D(1:len_C_Al-1))./2;

	d1 = [D_minus(1:end-1)*dtdx2, 0];
	d2 = [1, 1 - (D_plus(2:end) + D_minus(1:end-1))*dtdx2, 1];
	d3 = [0, D_plus(2:end)*dtdx2];
	Matrix = diag(d1, -1) + diag(d2) + diag(d3, +1);

	C_Al = Matrix*C_Al;
	C_Si = Matrix*C_Si;
	ind = find(j == checkTime); 
	if ind
		subplot(2, 1, 1);
		plot(1:len_C_Al, C_Al'/n_Al);
		hold on;

		subplot(2, 1, 2);
		plot( 1:len_C_Al, C_Si');
		hold on;
	end
end
C_Si'