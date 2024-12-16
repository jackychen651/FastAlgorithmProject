clear;clc;close all;
% Parameters
kappa = 2 * pi;           % Wave number (adjust as needed)
b = 0.5;                  % Scattering coefficient
R = 2;                  % Radius of the circular scatterer
L = 10;                   % Domain extends from -L/2 to L/2 in both x and y
Nx = 500;                 % Number of grid points in x
Ny = 500;                 % Number of grid points in y

% Grid spacing
dx = L / (Nx - 1);
dy = L / (Ny - 1);

% Generate grid
x = linspace(-L/2, L/2, Nx);
y = linspace(-L/2, L/2, Ny);
[X, Y] = meshgrid(x, y);

% Polar coordinates for scatterer
r = sqrt(X.^2 + Y.^2);
theta = atan2(Y, X);

% Incident field
u_in = exp(1i * kappa * X);

% Define scatterer mask
scatterer = r <= R;

% Define PML parameters
pml_width = 20;  % Number of grid points in PML
sigma_max = 25;  % Maximum damping coefficient in PML

% Create damping profiles for PML
sigma_x = zeros(Nx, 1);
sigma_y = zeros(Ny, 1);

% Damping in x-direction
for i = 1:Nx
    if i <= pml_width
        sigma_x(i) = sigma_max * ((pml_width - i + 1)/pml_width)^2;
    elseif i > Nx - pml_width
        sigma_x(i) = sigma_max * ((i - (Nx - pml_width))/pml_width)^2;
    else
        sigma_x(i) = 0;
    end
end

% Damping in y-direction
for j = 1:Ny
    if j <= pml_width
        sigma_y(j) = sigma_max * ((pml_width - j + 1)/pml_width)^2;
    elseif j > Ny - pml_width
        sigma_y(j) = sigma_max * ((j - (Ny - pml_width))/pml_width)^2;
    else
        sigma_y(j) = 0;
    end
end

% Create 2D damping coefficients
[Sig_x, Sig_y] = meshgrid(sigma_x, sigma_y);

% Total damping
Sigma = Sig_x + Sig_y;

% Total number of grid points
N = Nx * Ny;

% Mapping from 2D to 1D index
index = @(i,j) (j-1)*Nx + i;

% Initialize sparse matrix and RHS vector
A = sparse(N, N);
f = -kappa^2 * b * u_in(:);

% Precompute coefficients
coeff_center = -2/dx^2 - 2/dy^2 + kappa^2 * (1 - b) + 1i * Sigma(:);
coeff_x = 1/dx^2;
coeff_y = 1/dy^2;

% Assemble matrix A
for j = 1:Ny
    for i = 1:Nx
        p = index(i,j);
        if scatterer(p)  % Inside the scatterer
            % Assume the scatterer imposes a Dirichlet condition (u_sc = 0)
            A(p, p) = 1;
            f(p) = 0;
        elseif i == 1 || i == Nx || j == 1 || j == Ny
            % Boundary condition with PML: Implement absorbing boundary
            A(p, p) = 1 + 1i * Sigma(p)*dx^2;  % Adjusted for PML
            f(p) = 0;
        else
            % Interior points
            A(p,p) = coeff_center(p);
            A(p, index(i+1,j)) = coeff_x;
            A(p, index(i-1,j)) = coeff_x;
            A(p, index(i,j+1)) = coeff_y;
            A(p, index(i,j-1)) = coeff_y;
        end
    end
end
% Solve the system
disp('Solving the linear system...');
u_scattered = A \ f;
disp('Solution completed.');

% Reshape the scattered field to 2D grid
u_scattered = reshape(u_scattered, Nx, Ny)';
% Total field
u_total = u_in + u_scattered;

% Plot magnitude of scattered field
figure;
imagesc(x, y, abs(u_scattered));
colorbar;
xlabel('x');
ylabel('y');
title('Magnitude of Scattered Field |u^{sc}|');
axis equal tight;
set(gca, 'YDir', 'normal');  % Correct y-axis direction

% Plot phase of scattered field
figure;
imagesc(x, y, angle(u_scattered));
colorbar;
xlabel('x');
ylabel('y');
title('Phase of Scattered Field arg(u^{sc})');
axis equal tight;
set(gca, 'YDir', 'normal');

% Plot total field magnitude
figure;
imagesc(x, y, abs(u_total));
colorbar;
xlabel('x');
ylabel('y');
title('Magnitude of Total Field |u|');
axis equal tight;
set(gca, 'YDir', 'normal');