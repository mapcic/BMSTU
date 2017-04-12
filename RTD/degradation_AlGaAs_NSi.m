clear; clc;

e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); me = 9.11*1e-31;
nm = 1e-9; 
hbar = 1.054*1e-34; k_B = 1.38e-23;

T = 300;
kT = T*k_B;
Time = 10*365*24;

% % Count layers
% a = 8*6;
% b = 4*6;
% c = 6*6;
% Count layers
a = 8;
b = 4;
c = 6;

% Atoms and addop count
n_Atoms = 4.42*1e28;
n_Al = n_Atoms/2;
purity = 1e-8;

dx = 0.56*nm; 
dt = 1;
dtdx2 = dt*60*60/dx^2; 

checkTime = [0, 0.2, 0.5, 0.8, 1:2:10]*365*24;
grid_n_Al = [zeros(1, a), n_Al*ones(1, b), zeros(1, c), n_Al*ones(1, b), zeros(1, a)];
len_grid = length(grid_n_Al);

for i = 0 : dt : Time
	[grid_Ec, grid_Eg, grid_me_eff, grid_mp_eff] = getBandPropAlGaAs(grid_n_Al);
	
	Nc = 2*(me*grid_me_eff*kT/pi/hbar^2/2).^(3/2);
	Nv = 2*(me*grid_mp_eff*kT/pi/hbar^2/2).^(3/2);

	ni = sqrt(Nc.*Nv).*exp(-grid_Eg./(2*kT*JtoEv));	

	grid_n_conduct = 0.5*n_Atoms*purity*( 2 + 0.5*(2*ni./(n_Atoms*purity)).^2 );

	D = 0.2*exp(-3.5/(kT*JtoEv))*1e-4*(1e14*1e6./ni).^3;

	D_plus = (D(1:len_grid-1) + D(2:len_grid))./2;
	D_minus = (D(2:len_grid) + D(1:len_grid-1))./2;
	
	d1 = D_plus*dtdx2;
	d2 = [1 - D_plus(1)*dtdx2, 1 - (D_plus(2:end) + D_minus(1:end-1))*dtdx2, 1 - D_minus(end)*dtdx2];
	d3 = D_minus*dtdx2;

	Matrix = diag(d1, -1) + diag(d2) + diag(d3, +1);
	
	grid_n_Al = (Matrix*grid_n_Al')';

	if find(i == checkTime)
		plot(1:1:a+b+c+b+a, grid_n_Al'./n_Al);
		hold on;
	end
end