clear; clc;

e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); 
me = 9.11*1e-31; nm = 1e-9;

% Fermi Energy
EFermi = 1.51*1e-20; % J

T = 800; % K

checkTime = [0, 20]; %years

% monolayer in heterostructure
hsSize = [15, 6, 4, 6, 15];
% hsSize = [2 2 1 2 2];
% grid Step
% hsStep = [1, 1, 0.01, 1, 1];
hsStep = [1, 1, 0.5, 1, 1];
% Concentration of Aliminum
hsAl = [0, 0.44, 0, 0.44, 0];

% atoms' radius
dx = 0.56; %nm

% Applyied voltage
dU = 0:0.01:0.5;

% Get profile Ec
[grids_Ec, grids_meff, grids_C_Al] = getDODS_AlGaAs(...
	dx*nm,...
	checkTime, ...
	T, ...
	hsSize, ...
	hsStep,...
	hsAl...
);

dxGrid = [];
xGrid = [];
C_Al = [];

for ind = 1 : length(hsSize)
	% C_Al = [ C_Al, hsAl(ind)*ones( size(1 : hsStep(ind) : hsSize(ind)) ) ];
	tmp = dx*hsStep(ind)*ones( 1, hsSize(ind)/hsStep(ind));
	C_Al = [C_Al, hsAl(ind)*tmp];
	dxGrid = [dxGrid, tmp];
end
% size(dxGrid)
xGrid = cumsum(dxGrid);
% dxGrid = dxGrid(1:end-1);
plot(xGrid, grids_C_Al);

% % get J from V
% for j = 1 : length(checkTime)
% 	% J(j, :) = getJ(dx*nm, ...
% 	% 	grids_meff(j, :)*me, ...
% 	% 	grids_Ec(j, :)*eVtoJ, ...
% 	% 	dU*eVtoJ, ...
% 	% 	EFermi...
% 	% );
% 	Trans(j, :) = getTransperent(dx*nm, ...
% 		grids_meff(j, :)*me, ...
% 		grids_Ec(j, :)*eVtoJ...
% 	);
% end

% % showResult();

% figure('units', 'normalized', 'outerposition', [0 0 1 1]);
% Axes = {
% 	subplot(1, 2, 1);
% 	subplot(1, 2, 2);
% };
% plot(Axes{1}, (0:a+b+c+b+a-1)*dx, grids_Ec);
% plot(Axes{2}, log(Trans)', (ones(length(checkTime), 1)*linspace(0, max(max(grids_Ec)), 10000))' );
% plot(Axes{2}, log(Trans)', (ones( length(checkTime), 1)*linspace(0, max(max(grids_Ec))))' );