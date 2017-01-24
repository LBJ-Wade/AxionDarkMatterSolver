close all
clear all
clc
% Philip Mocz (2017), Harvard University
% Solve 3D Schrodinger + Poisson equation in fourier space in periodic box
% kick-drift-kick 2nd order pseudo-spectral unitary method
% i hbar d_t psi = [-hbar^2/(2m) nabla^2 + V + g|psi|^2] psi

% boxsize = 1
% hbar = 1
% <rho> = 1
% G = ?

N = 128;64;  % resolution
G = 10000;   % XXX
t = 0;
Tfinal = 10;
dt = 0.0001; % XXX may wish to implement adaptive, cfl based value
g = 0; -0.1;1;0;


xlin = linspace(0,1,N+1)';  % Note, x=0 and x=1 are the same point!
xlin = xlin(1:N);           % so crop it
dx = xlin(2) - xlin(1);

[x, y, z] = meshgrid(xlin, xlin, xlin);

% IC
sigma = 0.03;
rho = 0.95 + 0.01*exp(-((x-0.5).^2+(y-0.5).^2+(z-0.5).^2)/2/sigma^2)/(sigma^3*sqrt(2*pi)^3);
rho = rho + 0.01*exp(-((x-0.2).^2+(y-0.7).^2+(z-0.4).^2)/2/sigma^2)/(sigma^3*sqrt(2*pi)^3);
rho = rho + 0.01*exp(-((x-0.4).^2+(y-0.6).^2+(z-0.6).^2)/2/sigma^2)/(sigma^3*sqrt(2*pi)^3);
rho = rho + 0.01*exp(-((x-0.6).^2+(y-0.8).^2+(z-0.6).^2)/2/sigma^2)/(sigma^3*sqrt(2*pi)^3);
rho = rho + 0.01*exp(-((x-0.8).^2+(y-0.2).^2+(z-0.4).^2)/2/sigma^2)/(sigma^3*sqrt(2*pi)^3);
sum(rho(:))/N^3
psi = sqrt(rho);
clear x;
clear y;
clear z;

% fourier space variables
klin = 2*pi*(-N/2:N/2-1)';
[kx, ky, kz] = meshgrid(klin, klin, klin);
kSq = kx.^2 + ky.^2 + kz.^2;
clear kx;
clear ky;
clear kz;

% Poisson solver
Vhat = -fftshift(fftn(4*pi*G*(rho-1))) ./ ( kSq  + (kSq==0));
V = ifftn(ifftshift(Vhat));
V = V - max(V(:));


fh = figure(1);
set(fh,'position',[0,0,500,500]);

%% time evolution
while t < Tfinal
    
    % potential - (1/2) kick
    psi = exp(-1.i*dt/2*(V + g*rho)).*psi;
    
    % kinetic - drift
    psihat = fftshift(fftn(psi));
    psihat = exp(dt * (-1.i*kSq/2)) .*psihat;
    psi = ifftn(ifftshift(psihat));
    
    rho = abs(psi).^2;
    
    % potential - (1/2) kick
    Vhat = -fftshift(fftn(4*pi*G*(rho-1))) ./ ( kSq  + (kSq==0));
    V = ifftn(ifftshift(Vhat));
    V = V - max(V(:));
    psi = exp(-1.i*dt/2*(V + g*rho)).*psi;
    
    t = t + dt
    
    % plot
    figure(1);
    imagesc(log10(rho(:,:,end/2+1)))
    caxis([-1 3])
    axis off
    pause(0.001);
end
