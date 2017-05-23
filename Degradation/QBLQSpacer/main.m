clear; clc;

e = 1.6e-19; eVtoJ = e; JtoEv = e^(-1); 
me = 9.11*1e-31; nm = 1e-9;

% Do it smooth
dis = 1;
% atoms' radius
dx = 0.56; %nm
dx = dx/dis;
 
% Count layers
% a = [3, 7, 10, 15]; % monolayers
% b = 6;
% c = 6;

a = [3, 7, 10, 15]; % monolayers
b = 5;
c = 6;


a = a*dis;
b = b*dis;
c = c*dis;

sizeHS = a + b + c + b + a;

% Fermi Energy
EFermi = 1.51*1e-20; % J

% Applyied voltage
dU = 0:0.1:0.6;

% Ec
Ec1 = [zeros(1, a(1)), 0.5*ones(1, b), zeros(1, c), 0.5*ones(1, b), zeros(1, a(1))];
Ec2 = [zeros(1, a(2)), 0.5*ones(1, b), zeros(1, c), 0.5*ones(1, b), zeros(1, a(2))];
Ec3 = [zeros(1, a(3)), 0.5*ones(1, b), zeros(1, c), 0.5*ones(1, b), zeros(1, a(3))];
Ec4 = [zeros(1, a(4)), 0.5*ones(1, b), zeros(1, c), 0.5*ones(1, b), zeros(1, a(4))];

% meff
meff1 = [0.067*ones(1, a(1)), 0.15*ones(1, b), 0.067*ones(1, c), 0.15*ones(1, b), 0.067*ones(1, a(1))];
meff2 = [0.067*ones(1, a(2)), 0.15*ones(1, b), 0.067*ones(1, c), 0.15*ones(1, b), 0.067*ones(1, a(2))];
meff3 = [0.067*ones(1, a(3)), 0.15*ones(1, b), 0.067*ones(1, c), 0.15*ones(1, b), 0.067*ones(1, a(3))];
meff4 = [0.067*ones(1, a(4)), 0.15*ones(1, b), 0.067*ones(1, c), 0.15*ones(1, b), 0.067*ones(1, a(4))];

reserve = 15;

% numPoint = 5000;

% Tr1 = getTransperent(...
% 	dx*nm, ...
% 	meff1*me, ...
% 	Ec1*eVtoJ,...
% 	numPoint...
% );
tic
J1 = getJ(dx*nm, ...
	meff1*me, ...
	Ec1*eVtoJ, ...
	dU*eVtoJ, ...
	EFermi, ...
	reserve, ...
	a(1),...
	b,...
	c...
);
% Tr2 = getTransperent(...
% 	dx*nm, ...
% 	meff2*me, ...
% 	Ec2*eVtoJ,...
% 	numPoint...
% );

% J2 = getJ(dx*nm, ...
% 	meff2*me, ...
% 	Ec2*eVtoJ, ...
% 	dU*eVtoJ, ...
% 	EFermi, ...
% 	reserve, ...
% 	a(2),...
% 	b,...
% 	c...
% );

% % Tr3 = getTransperent(...
% % 	dx*nm, ...
% % 	meff3*me, ...
% % 	Ec3*eVtoJ,...
% % 	numPoint...
% % );

% J3 = getJ(dx*nm, ...
% 	meff3*me, ...
% 	Ec3*eVtoJ, ...
% 	dU*eVtoJ, ...
% 	EFermi, ...
% 	reserve, ...
% 	a(3),...
% 	b,...
% 	c...
% );


% % Tr4 = getTransperent(...
% % 	dx*nm, ...
% % 	meff4*me, ...
% % 	Ec4*eVtoJ,...
% % 	numPoint...
% % );

% J4 = getJ(dx*nm, ...
% 	meff4*me, ...
% 	Ec4*eVtoJ, ...
% 	dU*eVtoJ, ...
% 	EFermi, ...
% 	reserve, ...
% 	a(4),...
% 	b,...
% 	c...
% );
toc

% showResult(dx*nm, sizeHS, Ec1, Ec2, Ec3, Ec4, J1, J2, J3, J4, dU, Tr1, Tr2, Tr3, Tr4, a);
showResult(dx*nm, sizeHS, Ec1, Ec2, Ec3, Ec4, J1, J2, J3, J4, dU, a);
