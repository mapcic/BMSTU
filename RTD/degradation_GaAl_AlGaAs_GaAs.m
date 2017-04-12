clear; clc;

e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); me = 9.1*1e-31;
nm = 1e-9; 
hbar = 1.054*1e-34; k_B = 1.38e-23;

T = 800;
kT = T*k_B;
Time = 25*365*24;

n_Atoms = 4.42*1e28;
n_Al = n_Atoms/2;

dx = 0.56*nm; 
dt = 1;
dtdx2 = dt*60*60/dx^2; 

D = 0.2*exp(-3.5/(kT*JtoEv))*1e-4;

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

plot(1:1:a+b+c+b+a, grid_n./n_Al);
hold on;

checkTime = [1, 5, 10, 15, 20, 25]*365*24;

grid_n = grid_n';
for i = 0 : dt : Time
	grid_n = Matrix*grid_n;
	if find(i == checkTime)
		plot(1:1:a+b+c+b+a, grid_n'./n_Al);
	end
end