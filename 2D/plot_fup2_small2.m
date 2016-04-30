clear all; close all; clc;

N = 2^8;
h = 4/(N-1);
x = -2:h:2;
y1 = fup2_small(x);

% plot
figure('color','w')
hold on
plot(x,y1);
LW = 'LineWidth';
lw = 2;
px = [-1 0 1];
plot(px,fup2_small(px),'ks',LW,lw)
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
info = ['$\mathrm{fup}_2(x) = \frac{36}{16}(' ...
    '\mathrm{up}(\frac{x}{4}-\frac{1}{2}) - 2\mathrm{up}(\frac{x}{4}-\frac{3}{4}) +' ...
    '2\mathrm{up}(\frac{x}{4}-1) - ' ...
    '2\mathrm{up}(\frac{x}{4}-\frac{5}{4})$'];
title(info,'interpreter','latex','fontsize',12)

figure('color','w')
% hold on
% plot(x,y);
d1y = diff(y1) / h;
d2y = diff(d1y) / h;
plot(x,y1, x(2:end),d1y, x(3:end),d2y);
title('$\mathrm{fup}_2, \mathrm{fup}''_2, \mathrm{fup}''''_2$','interpreter','latex','fontsize',12)