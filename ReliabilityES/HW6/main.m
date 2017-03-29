clear; clc;

syms t;

f1 = 1 + 0*t;
f2 = 0 + 2*t;
f3 = 2 - 2*t;

F1 = 1 - int(f1, t);
F2 = 1 - int(f2, t);
F3 = 1 - int(f3, t);

lmbd1 = f1./F1;
lmbd2 = f2./F2;
lmbd3 = f3./F3;

calc = @(vec, func) double(subs(func, t, vec));

dt = 1e-2;
T = dt*1e-1 : dt : 1-dt*1e-1;

lmbd1 = calc(T, lmbd1);
lmbd2 = calc(T, lmbd2);
lmbd3 = calc(T, lmbd3);
lmbd0 = lmbd1 + lmbd2 + lmbd3;

windows = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
Axes = {
    subplot(2, 2, 1)
    subplot(2, 2, 2)
    subplot(2, 2, 3)
    subplot(2, 2, 4)
};

plot(Axes{1}, T, lmbd1);
plotFormat( Axes{1}, '$\lambda_{1}$', '$t,\,h$', '$\lambda,\,h^{-1}$', {}, [0.7, 1], [0, 100], 16);

plot(Axes{2}, T, lmbd2);
plotFormat( Axes{2}, '$\lambda_{2}$', '$t,\,h$', '$\lambda,\,h^{-1}$', {}, [0.7, 1], [0, 100], 16);

plot(Axes{3}, T, lmbd3);
plotFormat( Axes{3}, '$\lambda_{3}$', '$t,\,h$', '$\lambda,\,h^{-1}$', {}, [0.7, 1], [0, 100], 16);

plot(Axes{4}, T, lmbd0);
plotFormat( Axes{4}, '$\lambda_{sys}$', '$t,\,h$', '$\lambda,\,h^{-1}$', {}, [0.7, 1], [0, 100], 16);