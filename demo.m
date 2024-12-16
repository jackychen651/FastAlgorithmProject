clear;clc;close all;
L=12;
h = 2.0^(-L);
omega = 25.0;

a = 1;
b = 1;

x = -a/2:h:a/2-h;
y = -b/2:h:b/2-h;

n = length(x); m = length(y);

[X,Y] = meshgrid(x,y);
nu = @(x,y) 1.5*exp(-160*(x.^2 + y.^2));

M = nu(X,Y);



tic
LS = LippmannSchwinger(x,y,omega,nu,a);

u_inc = exp(omega*1i*X);
rhsDual = -omega^2*nu(X,Y).*u_inc ;

sigma = LS\rhsDual(:);

u = LS.apply_Green(sigma);


figure(1); clf();
imagesc(M);
colorbar;
filename = sprintf('Figure1_L=%d.pdf', L);
set(gcf, 'PaperUnits', 'centimeters');
paperWidth = 15;
paperHeight = 10;
set(gcf, 'PaperSize', [paperWidth paperHeight]);
set(gcf, 'PaperPosition', [0 0 paperWidth paperHeight]);
saveas(gcf, filename);

figure(2); clf();
imagesc(real(reshape(u, n, m)));
colorbar;
filename = sprintf('Figure2_L=%d.pdf', L);
set(gcf, 'PaperUnits', 'centimeters');
paperWidth = 15;
paperHeight = 10;
set(gcf, 'PaperSize', [paperWidth paperHeight]);
set(gcf, 'PaperPosition', [0 0 paperWidth paperHeight]);
saveas(gcf, filename);

figure(3); clf();
imagesc(real(reshape(u + u_inc, n, m)));
colorbar;
filename = sprintf('Figure3_L=%d.pdf', L);
set(gcf, 'PaperUnits', 'centimeters');
paperWidth = 15;
paperHeight = 10;
set(gcf, 'PaperSize', [paperWidth paperHeight]);
set(gcf, 'PaperPosition', [0 0 paperWidth paperHeight]);
saveas(gcf, filename);


colorbar
toc