clear; clc;

e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); 
me = 9.11*1e-31; nm = 1e-9;

% Do it smooth
dis = 1;
% atoms' radius
dx = 0.56; %nm
dx = dx/dis;

% Count layers
a = 4; % monolayers
b = 3;
c = 5;

a = a*dis;
b = b*dis;
c = c*dis;

sizeBHS = a + b + a;
sizeRTHS = a + b + c + b + a;

% Fermi Energy
EFermi = 1.51*1e-20; % J

% Applyied voltage
dU = 0:0.001:0.5;

% Ec
EcBHS = [zeros(1, a), 1*ones(1, b), zeros(1, a)];
EcRTHS = [zeros(1, a), 1*ones(1, b), zeros(1, c), 1*ones(1, b), zeros(1, a)];

% meff
meffBHS = [0.067*ones(1, a), 0.15*ones(1, b), 0.067*ones(1, a)];
meffRTHS = [0.067*ones(1, a), 0.15*ones(1, b), 0.067*ones(1, c), 0.15*ones(1, b), 0.067*ones(1, a)];

numPoint = 1000;

TrBHS = getTransperent(...
	dx*nm, ...
	meffBHS*me, ...
	EcBHS*eVtoJ,...
	numPoint...
);

TrRTHS = getTransperent(...
	dx*nm, ...
	meffRTHS*me, ...
	EcRTHS*eVtoJ,...
	numPoint...
);

JBHS = getJ(dx*nm, ...
	meffBHS*me, ...
	EcBHS*eVtoJ, ...
	dU*eVtoJ, ...
	EFermi...
);

JRTHS = getJ(dx*nm, ...
	meffRTHS*me, ...
	EcRTHS*eVtoJ, ...
	dU*eVtoJ, ...
	EFermi...
);

showResult(dx, EcRTHS, EcBHS, JBHS, JRTHS, TrBHS, TrRTHS, dU);