clear; clc;

e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1);
nm = 1e-9; me = 9.1*1e-31;
hbar = 1.054*1e-34; k_B = 1.38e-23;

T = 300+273.5;
kT = T*k_B; % J

% Count layers
a = 8; % monolayers
b = 5;
c = 6;

Time = 28*60*60; % to hours

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
dtdx2 = dt/dx^2; % s/m^2

D_Al = 0.2*exp(-3.5/(kT*JtoEv))*1e-4; % m^2/s
D_Si = 0.17*exp(-3.5/(kT*JtoEv))*1e-4; % m^2/s
% D = 2.9*1e8*exp(-6/(kT*JtoEv))*1e-4;

n_Si = 2*1e17*1e6;
C_Si = [n_Si, ni*ones(1, a), ni*ones(1, b), ni*ones(1, c), ni*ones(1, b), ni*ones(1, a), n_Si];

x_Al = [0, zeros(1, a), ones(1, b), zeros(1, c), ones(1, b), zeros(1, a), 0];
C_Al = x_Al*n_Al;
len_C_Al = length(C_Al);

d1 = D_Al*dtdx2*ones(1, len_C_Al-1);
d2 = [ 1 - D_Al*dtdx2, 1 - 2*D_Al*dtdx2*ones(1, len_C_Al-2), 1 - D_Al*dtdx2 ];
d3 = D_Al*dtdx2*ones(1, len_C_Al-1);
Matrix_Al = diag(d1, -1) + diag(d2) + diag(d3, +1);

d1 = [D_Si*dtdx2*ones(1, len_C_Al-2), 0];
d2 = [ 1, 1 - 2*D_Si*dtdx2*ones(1, len_C_Al-2), 1 ];
d3 = [0, D_Si*dtdx2*ones(1, len_C_Al-2)];
Matrix_Si = diag(d1, -1) + diag(d2) + diag(d3, +1);

checkTime = (0:28)*60*60;

C_Al = C_Al';
C_Si = C_Si';
for j = 0 : dt : Time
	D_Si = 0.2*exp(-3.5/(kT*JtoEv))*(C_Si'./ni).^3*1e-4;

	D_plus = (D_Si(1:len_C_Al-1) + D_Si(2:len_C_Al))./2;
	D_minus = (D_Si(2:len_C_Al) + D_Si(1:len_C_Al-1))./2;

	d1 = [D_minus(1:end-1)*dtdx2, 0];
	d2 = [1, 1 - (D_plus(2:end) + D_minus(1:end-1))*dtdx2, 1];
	d3 = [0, D_plus(2:end)*dtdx2];
	Matrix_Si = diag(d1, -1) + diag(d2) + diag(d3, +1);

	D_Al = 0.2*exp(-3.5/(kT*JtoEv))*(C_Si'./ni).^3*1e-4;

	D_plus = (D_Al(1:len_C_Al-1) + D_Al(2:len_C_Al))./2;
	D_minus = (D_Al(2:len_C_Al) + D_Al(1:len_C_Al-1))./2;

	d1 = [D_minus(1:end-1)*dtdx2, 0];
	d2 = [1, 1 - (D_plus(2:end) + D_minus(1:end-1))*dtdx2, 1];
	d3 = [0, D_plus(2:end)*dtdx2];
	Matrix_Al = diag(d1, -1) + diag(d2) + diag(d3, +1);

	C_Al = Matrix_Al*C_Al;
	C_Si = Matrix_Si*C_Si;
	ind = find(j == checkTime); 
	if ind
		% plot( 1:len_C_Al, C_Al');
		% hold on;
		subplot(2, 1, 1);
		plot(1:len_C_Al, C_Al');
		hold on;

		subplot(2, 1, 2);
		plot( 1:len_C_Al, C_Si'./ni);
		hold on;
		% [Ec(ind, :), meff(ind, :)] = getEcAlGaAs(C_Al');
	end
end