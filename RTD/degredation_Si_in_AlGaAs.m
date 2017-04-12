clear; clc;

e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); me = 9.11*1e-31;
nm = 1e-9; 
hbar = 1.054*1e-34; k_B = 1.38e-23;

T = 400;
kT = T*k_B;
Time = 5*24*60 + 23;

% % Count layers. Smooth
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
purity = 1e-9;

dx = 0.56*nm; 
dt = 1;
dtdx2 = dt*60/dx^2; 

n_Si = 2*1e18*1e6;
n_added = n_Atoms*purity;

% dtdx2 = dt*60*60/dx^2;

grid_n_Si = [n_Si, n_added*ones(1, a), n_added*ones(1, b), n_added*ones(1, c), n_added*ones(1, b), n_added*ones(1, a), n_Si];

grid_U = [0, zeros(1, a), ones(1, b), zeros(1, c), ones(1, b), zeros(1, a), 0];

len_grid = length(grid_n_Si);

grid_n_Al = [0, zeros(1, a), n_Al*ones(1, b), zeros(1, c), n_Al*ones(1, b), zeros(1, a), 0];
[grid_Ec, grid_Eg, grid_me_eff, grid_mp_eff] = getBandPropAlGaAs(grid_n_Al);
Nc = 2*(me*grid_me_eff*kT/pi/hbar^2/2).^(3/2);
Nv = 2*(me*grid_mp_eff*kT/pi/hbar^2/2).^(3/2);

ni = sqrt(Nc.*Nv).*exp(-grid_Eg./(2*kT*JtoEv));	

checkTime = (0:1:25)*365*24;
for i = 0 : dt : Time
	grid_n_conduct = 0.5*grid_n_Si.*( 2 + 0.5*(2*ni./(grid_n_Si)).^2 );

	D = 0.2*exp(-3.5/(kT*JtoEv))*1e-4*(grid_n_conduct./ni).^3;

	D_plus = (D(1:len_grid-1) + D(2:len_grid))./2;
	D_minus = (D(2:len_grid) + D(1:len_grid-1))./2;
	
	d1 = [D_plus(1:end-1)*dtdx2, 0];
	d2 = [1, 1 - (D_plus(2:end) + D_minus(1:end-1))*dtdx2, 1];
	d3 = [0, D_minus(2:end)*dtdx2];

	Matrix = diag(d1, -1) + diag(d2) + diag(d3, +1);
	
	grid_n_Si = (Matrix*grid_n_Si')';
	% if find(i == checkTime)
	% 	% i/365/24/60/60
	% 	plot(0:1:a+b+c+b+a+1, grid_n_Si'./n_Si);
	% 	hold on;
	% 	if ( isnan(grid_n_Si ))
	% 		10101
	% 	end
	% end
end

grid_n_Si

% plot(0:1:a+b+c+b+a+1, grid_n_Si'./n_Si, '--');
% hold on;
% plot(0:1:a+b+c+b+a+1, grid_U);