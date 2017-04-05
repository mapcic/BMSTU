clear; clc;

e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); me = 9.1*1e-31;
nm = 1e-9; 
hbar = 1.054*1e-34; k_B = 1.38e-23;

T = 350;
kT = T*k_B;
Time = 100*365*24*60;

Eg = 1.519 - 5.405*1e-4*T^2/(T+204);

me_eff = 0.067;
mp_eff = 0.51;

Nc = 2*(me*me_eff*kT/pi/hbar^2/2)^(3/2);
Nv = 2*(me*mp_eff*kT/pi/hbar^2/2)^(3/2);

ni = sqrt(Nc*Nv)*exp(-Eg/(2*kT*JtoEv));

n_Al = 4.42*1e28;
purity = 1e-10;

dx = 0.56*nm; 
dt = 1;
dtdx2 = dt*60/dx^2; 

n_conduct = 0.5*n_Al*purity*( 2 + 0.5*(2*ni/(n_Al*purity))^2 );

% D = 0.2*exp(-3.5/(kT*JtoEv))*1e-4;
D = 0.2*exp(-3.5/(kT*JtoEv))*1e-4*(n_conduct/ni)^3;

% Count layers
a = 8;
b = 4;
c = 6;

grid_n = [zeros(1, a), n_Al*ones(1, b), zeros(1, c), n_Al*ones(1, b), zeros(1, a)];
len_grid_n = length(grid_n);

d1 = [D*dtdx2*ones(1, len_grid_n-1)];
d2 = [1 - D*dtdx2, 1 - 2*D*dtdx2*ones(1, len_grid_n-2), 1 - D*dtdx2 ];
d3 = [D*dtdx2*ones(1, len_grid_n-1)];

Matrix = diag(d1, -1) + diag(d2) + diag(d3, +1);

plot(1:1:a+b+c+b+a, grid_n);
hold on;

grid_n = grid_n';
for i = 0 : dt : Time
	grid_n = Matrix*grid_n;
end

plot(1:1:a+b+c+b+a, grid_n');