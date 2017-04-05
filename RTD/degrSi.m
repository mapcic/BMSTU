clear; clc;

e = 1.6e-19;
eVtoJ = e;
JtoEv = e^(-1); 
nm = 1e-9;
me = 9.10938356*1e-31;
hbar = 1.054*1e-34;
k_B = 1.38e-23;
T = 400;
kT = T*k_B;

dx = 0.56*nm;
dt = 1;
Time = 15*365*24;

Eg = 1.519 - 5.405*1e-4*T^2/(T+204);

me_eff = 0.067;
mp_eff = 0.51;

Nc = 2*(me*me_eff*kT/pi/hbar^2/2)^(3/2);
Nv = 2*(me*mp_eff*kT/pi/hbar^2/2)^(3/2);

ni = sqrt(Nc*Nv)*exp(-Eg/(2*kT*JtoEv));

% Count layers
a = 8;
b = 4;
c = 6;


n_Si = 2*1e17*1e6;
grid_n = [n_Si, ni*ones(1, a), ni*ones(1, b), ni*ones(1, c), ni*ones(1, b), ni*ones(1, a), n_Si];
len_grid_n = length(grid_n);

grid_U = [0, zeros(1, a), n_Si*ones(1, b), zeros(1, c), n_Si*ones(1, b), zeros(1, a), 0];

plot(0:1:a+b+c+b+a+1, grid_n);
hold on;

for i = 0 : dt : Time
	dtdx2 = dt*60*60/dx^2;

	D = 0.2*exp(-3.5/(kT*JtoEv))*(grid_n./ni).^3*1e-4;
	D_plus = (D(1:len_grid_n-1) + D(2:len_grid_n))./2;
	D_minus = (D(2:len_grid_n) + D(1:len_grid_n-1))./2;
	
	d1 = [D_minus(1:end-1)*dtdx2, 0];
	d2 = [1, 1 - (D_plus(2:end) + D_minus(1:end-1))*dtdx2, 1];
	d3 = [0, D_plus(2:end)*dtdx2];

	Matrix = diag(d1, -1) + diag(d2) + diag(d3, +1);
	grid_n = (Matrix*grid_n')';

end

plot(0:1:a+b+c+b+a+1, grid_n', '--');
plot(0:1:a+b+c+b+a+1, grid_U);
