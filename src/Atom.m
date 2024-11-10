% Constants
hbar = 1.055e-34; % Planck's constant (J.s)
e = 1.602e-19; % Elementary charge (C)
epsilon0 = 8.854e-12; % Vacuum permittivity (F/m)
m_e = 9.109e-31; % Electron mass (kg)

mu = m_e / 2; % Reduced mass
r_min = 1e-10; % Minimum radius (m)
r_max = 20e-10; % Maximum radius (m)
N_r = 500; % Number of radial points
r = linspace(r_min, r_max, N_r); % Radial grid

V = -e^2 ./ (4 * pi * epsilon0 * r); % Coulomb potential

h = r(2) - r(1); % Step size
T = -hbar^2 / (2 * mu * h^2);
H_kinetic = diag(-2 * ones(N_r, 1)) + diag(ones(N_r-1, 1), 1) + diag(ones(N_r-1, 1), -1);
H_kinetic = T * H_kinetic;

H_potential = diag(V);
H = H_kinetic + H_potential;

[eigenvectors, eigenvalues] = eig(H);
energies = diag(eigenvalues);

[energies, idx] = sort(energies);
wavefunctions = eigenvectors(:, idx);

psi_ground = wavefunctions(:, 1);
psi_ground_normalized = psi_ground / sqrt(trapz(r, abs(psi_ground).^2));

probability_density = abs(psi_ground_normalized).^2;

% Plotting the probability density
figure;
plot(r, probability_density, 'LineWidth', 2);
xlabel('Radius (m)');
ylabel('Probability Density');
title('Hydrogen Atom Ground State Orbital (1s)');
grid on;

fprintf('Ground state energy: %.4e J\n', energies(1));
fprintf('Ground state energy: %.4f eV\n', energies(1) / e);

% Visualization of the orbital shape
% Create a meshgrid for spherical coordinates
theta = linspace(0, pi, 50);
phi = linspace(0, 2*pi, 50);
[Theta, Phi] = meshgrid(theta, phi);

% Calculate the radial part for the ground state (1s orbital)
R_1s = (1/sqrt(pi)) * exp(-r_min); % Normalized radial wave function for 1s

% Convert spherical coordinates to Cartesian
X = R_1s * sin(Theta) .* cos(Phi);
Y = R_1s * sin(Theta) .* sin(Phi);
Z = R_1s * cos(Theta);

% Plotting the orbital shape
figure;
surf(X, Y, Z, 'EdgeColor', 'none');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Plot of Hydrogen Atom Ground State Orbital (1s)');
axis equal;
view(30, 30); % Adjust view angle for better visualization
colorbar;
