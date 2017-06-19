clear; clc;

e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); 
me = 9.11*1e-31; nm = 1e-9;
hbar = 1.054*1e-34; k_B = 1.38e-23;

% Do it smooth
dis = 1;
% atoms' radius
dx = 0.56/2; %nm
dx = dx/dis;

% Count layers
a = 10*2; % monolayers
b = 8*2;
c = 5*2;
r = 30*2;

a = a*dis;
b = b*dis;
c = c*dis;

sizeBHS = a + b + a;
sizeRTHS = a + b + c + b + a;

% Fermi Energy
EFermi = 1.51*1e-20; % J

% Applyied voltage
% dU = 0:0.001:0.5;
dU = [0, 0.25, 0.5];

% Ec
EcRTHS = [zeros(1, a), 1*ones(1, b), zeros(1, c), 1*ones(1, b), zeros(1, a)];

% meff
meffRTHS = [0.067*ones(1, a), 0.15*ones(1, b), 0.067*ones(1, c), 0.15*ones(1, b), 0.067*ones(1, a)];


ni_GaAs = 1e12;
ni_AlAs = 1e12;
Nd = 1e24;
Ni = [Nd*ones(1, r), ni_GaAs*ones(1, a), ni_AlAs*ones(1, b), ni_GaAs*ones(1, c), ni_AlAs*ones(1, b), ni_GaAs*ones(1, a), Nd*ones(1, r)];

eps = 13.18 - 3.12*[zeros(1, r), zeros(1, a), ones(1, b), zeros(1, c), ones(1, b), zeros(1, a), zeros(1, r)];

boundL = r;
boundR = r + length(EcRTHS);

for j = 1:length(dU)
	[V1(j, :), n1(j, :)] = getConcentrationElectrons(...
		0.001,...
		[EcRTHS(1)*ones(1, r), EcRTHS, EcRTHS(end)*ones(1, r)]*eVtoJ,...
		[meffRTHS(1)*ones(1, r), meffRTHS, meffRTHS(end)*ones(1, r)]*me,...
		Ni,...
		eps,...
		dx*nm,...
		dU(j)*eVtoJ,...
		boundL,...
		boundR);
end


ni_GaAs = 1e12;
ni_AlAs = 1e7;
Nd = 1e24;
Ni = [Nd*ones(1, r), ni_GaAs*ones(1, a), ni_AlAs*ones(1, b), ni_GaAs*ones(1, c), ni_AlAs*ones(1, b), ni_GaAs*ones(1, a), Nd*ones(1, r)];

eps = 13.18 - 3.12*[zeros(1, r), zeros(1, a), ones(1, b), zeros(1, c), ones(1, b), zeros(1, a), zeros(1, r)];

boundL = r;
boundR = r + length(EcRTHS);

for j = 1:length(dU)
	[V2(j, :), n2(j, :)] = getConcentrationElectrons(...
		0.001,...
		[EcRTHS(1)*ones(1, r), EcRTHS, EcRTHS(end)*ones(1, r)]*eVtoJ,...
		[meffRTHS(1)*ones(1, r), meffRTHS, meffRTHS(end)*ones(1, r)]*me,...
		Ni,...
		eps,...
		dx*nm,...
		dU(j)*eVtoJ,...
		boundL,...
		boundR);
end

plot(1:length(EcRTHS), -V1(1, r+1:length(EcRTHS)+r) + EcRTHS, 1:length(EcRTHS), -V2(1, r+1:length(EcRTHS)+r) + EcRTHS);
% plot(1:length(n(1, :)), n(3, :));

% showResult(dx, EcRTHS, EcBHS, JBHS, JRTHS, TrBHS, TrRTHS, dU);