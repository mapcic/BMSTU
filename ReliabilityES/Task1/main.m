clear; clc;

% Include common function
addpath(strcat(fileparts(mfilename('fullpath')), '/../../common/'));

% Main part
N = 1000;  % N is amount of all irrecoverable products
n = [0 50 40 32 25 20 17 16 16 15 14 15 14 14 13 14 13 13 13 14 12 12 13 12 13 14 16 20 25 30 40]; % n is amount of breaked products at every moment

Nw = N - cumsum(n);
nLen = length(n);

dt = 100;
t = 100 : dt : 3000;
tLen = length(t);

f = n(2:nLen)./(dt*N);
F = Nw(2:nLen)./N;
lambda = 2.*n(2:nLen)./(dt.*( Nw(2:nLen) + Nw(1:nLen-1) ));

% Draw graphics
windows = figure('units', 'normalized', 'outerposition', [0 0 1 1]);

Axes = {
	subplot(1, 3, 1);
	subplot(1, 3, 2);
	subplot(1, 3, 3);
};

plotFormat(Axes{1}, 'Title', 'xTitle', 'yTitle', {}, [], []);
plot(Axes{1}, t, F);
bar(Axes{1}, t-dt/2, F);

plotFormat(Axes{2}, 'Title', 'xTitle', 'yTitle', {}, [], []);
plot(Axes{2}, t, f);
bar(Axes{2}, t-dt/2, f);

plotFormat(Axes{3}, 'Title', 'xTitle', 'yTitle', {}, [], []);
plot(Axes{3}, t, lambda);
bar(Axes{3}, t-dt/2, lambda);