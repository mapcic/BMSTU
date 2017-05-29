clear; clc;

e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); 
me = 9.11*1e-31; nm = 1e-9;

T = 800; % K

% checkTime = 0:25:25; %years

checkTime = [0]; %years

% atoms' radius
dx = 0.56; %nm
% Smooth
% dx = 0.3; %nm

% Count layers
% Active field
a = 10; % monolayers
b = 6;
c = 6;

sizeHS = a + b + c + b + a;

% Smooth
% a = a*2; % monolayers
% b = b*2;
% c = c*2;
% 
% sizeHS = a + b + c + b + a;

% Fermi Energy
EFermi = 1.51*1e-20; % J

% Applyied voltage
dU = 0:0.01:0.5;

% grid of Al conentration
grid_x_Al = [zeros(1, a), ...
	0.44*ones(1, b), ...
	zeros(1, c), ...
	0.44*ones(1, b), ...
	zeros(1, a)
];

% Get profile Ec
% [grids_Ec, grids_meff, grids_C_Al] = getDiffCloseAlGaAs( grid_x_Al, checkTime, dx*nm, T );
[grids_Ec, grids_meff, grids_C_Al] = getDiffOpenAlGaAs( grid_x_Al, checkTime, dx*nm, T );
% [grids_Ec, grids_meff, grids_C_Al] = getDiffCloseAlGaAsNd( grid_x_Al, checkTime, dx*nm, T, 5e15*1e6 );
% [grids_Ec, grids_meff, grids_C_Al] = getDiffOpenAlGaAsNd( grid_x_Al, checkTime, dx*nm, T, 5e15*1e6 );
% [grids_Ec, grids_meff, grids_C_Al, Six] = getDiffAlGaAs_Si( grid_x_Al, checkTime, dx*nm, T, 2*1e18*1e6 );

% get J from V
for j = 1 : length(checkTime)
	% J(j, :) = getJ(dx*nm, ...
	% 	grids_meff(j, :)*me, ...
	% 	grids_Ec(j, :)*eVtoJ, ...
	% 	dU*eVtoJ, ...
	% 	EFermi...
	% );
	Trans(j, :) = getTransperent(...
		dx*nm, ...
		grids_meff(j, :)*me, ...
		grids_Ec(j, :)*eVtoJ...
	);
end

% showResult();

figure('units', 'normalized', 'outerposition', [0 0 1 1]);
Axes = {
	subplot(1, 2, 1);
	subplot(1, 2, 2);
};
plot(Axes{1}, (0:a+b+c+b+a-1)*dx, grids_Ec);
plot(Axes{2}, log(Trans)', (ones(length(checkTime), 1)*linspace(0, max(max(grids_Ec)), 10000))' );
% plot(Axes{2}, log(Trans)', (ones( length(checkTime), 1)*linspace(0, max(max(grids_Ec))))' );