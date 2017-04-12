clear; clc;

e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); me = 9.11*1e-31;
nm = 1e-9;
hbar = 1.054*1e-34; k_B = 1.38e-23;

dx = 0.56*nm;
dt = 1;

T = 400;
kT = T*k_B;

Eg = 1.519 - 5.405*1e-4*T^2/(T+204);

me_eff = 0.067;
mp_eff = 0.51;

Nc = 2*(me*me_eff*kT/pi/hbar^2/2)^(3/2);
Nv = 2*(me*mp_eff*kT/pi/hbar^2/2)^(3/2);

ni = sqrt(Nc*Nv)*exp(-Eg/(2*kT*JtoEv));

C_Si = 2*1e17*1e6;

C_Al = 4.42*1e28;
purity = 1e-9;

C_addition = C_Al*purity;
% C_addition = ni;

Time = 10*365*24;
% Count layers
a = 8;
b = 4;
c = 6;

grid_C = [C_Si, C_addition*ones(1, a), C_addition*ones(1, b), C_addition*ones(1, c), C_addition*ones(1, b), C_addition*ones(1, a), C_Si];
len_grid_C = length(grid_C);

grid_U = [0, zeros(1, a), C_Si*ones(1, b), zeros(1, c), C_Si*ones(1, b), zeros(1, a), 0];

plot(0:1:a+b+c+b+a+1, len_grid_C);
hold on;

dtdx2 = dt*60*60/dx^2;
for i = 0 : dt : Time

	grid_n_coduct = 0.5*grid_C.*( 2 + 0.5*(2*ni./grid_C).^2 );

	D = 0.2*exp(-3.5/(kT*JtoEv))*(grid_n_coduct./ni).^3*1e-4;

	D_plus = (D(1:len_grid_C-1) + D(2:len_grid_C))./2;
	D_minus = (D(2:len_grid_C) + D(1:len_grid_C-1))./2;
	
	d1 = [D_minus(1:end-1)*dtdx2, 0];
	d2 = [1, 1 - (D_plus(2:end) + D_minus(1:end-1))*dtdx2, 1];
	d3 = [0, D_plus(2:end)*dtdx2];

	Matrix = diag(d1, -1) + diag(d2) + diag(d3, +1);
	
	grid_C = (Matrix*grid_C')';
end

plot(0:1:a+b+c+b+a+1, grid_C', '--');
plot(0:1:a+b+c+b+a+1, grid_U);