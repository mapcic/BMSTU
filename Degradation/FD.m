clear; clc;

windows = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
Axes = {
   subplot(1, 1, 1);
};
legends = {
  'I am empty...'  
};

mu = 100;
kT = [1 10 50 100];
E = 0:1:5*mu;

muMat = mu*ones(length(kT), length(E));
EMat = ones(length(kT), 1)*E;
kTMat = kT'*ones(1, length(E));

f_FD = 1./(1+exp((EMat-muMat)./(kTMat)));

kT = fliplr(kT);
for j = 1:length(kT)
   legends{j} = strcat('$kT = \frac{E_{F}}{', num2str(kT(j)), '}$'); 
end

plot(Axes{1}, E'/mu, f_FD', 'lineWidth', 2);
plotFormat(Axes{1}, 'Fermi-Dirac-Statistik, $f_{FD}(E)$', '$E/E_{F}$', '$\overline{n}$', legends, [], []);