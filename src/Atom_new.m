hbar = 1.055e-34;
e = 1.602e-19;
epsilon0 = 8.854e-12;
m_e = 9.109e-31;

atom = input('Enter the atom to visualize (H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca): ', 's');
orbital = input('Enter the orbital to visualize (1s, 2s, 2p, 3s, 3p, 3d, 4s): ', 's');

atomic_numbers = struct('H', 1, 'He', 2, 'Li', 3, 'Be', 4, 'B', 5, ...
                         'C', 6, 'N', 7, 'O', 8, 'F', 9, 'Ne', 10,...
                         'Na', 11, 'Mg', 12, 'Al', 13, 'Si', 14,...
                         'P', 15, 'S', 16, 'Cl', 17, 'Ar', 18,...
                         'K', 19, 'Ca', 20);

if isfield(atomic_numbers, atom)
    Z = atomic_numbers.(atom);
    mu = m_e;
else
    error('Invalid atom selected.');
end

r_min = 1e-10; 
r_max = 20e-10; 
N_r = 500; 
r = linspace(r_min, r_max, N_r); 

V = -Z * e^2 ./ (4 * pi * epsilon0 * r); 

h = r(2) - r(1); 
T = -hbar^2 / (2 * mu * h^2);
H_kinetic = diag(-2 * ones(N_r, 1)) + diag(ones(N_r-1, 1), 1) + diag(ones(N_r-1, 1), -1);
H_kinetic = T * H_kinetic;

H_potential = diag(V);
H = H_kinetic + H_potential;

[eigenvectors, eigenvalues] = eig(H);
energies = diag(eigenvalues);
[energies, idx] = sort(energies);
wavefunctions = eigenvectors(:, idx);

switch orbital
    case '1s'
        psi = wavefunctions(:, 1);
    case '2s'
        psi = wavefunctions(:, min(2,Z));
    case '2p'
        psi = wavefunctions(:, min(3,Z));
    case '3s'
        psi = wavefunctions(:, min(4,Z));
    case '3p'
        psi = wavefunctions(:, min(5,Z));
    case '3d'
        psi = wavefunctions(:, min(6,Z));
    case '4s'
        psi = wavefunctions(:, min(7,Z));
    otherwise
        error('Invalid orbital selected.');
end

psi_normalized = psi / sqrt(trapz(r, abs(psi).^2));
probability_density = abs(psi_normalized).^2;

figure;
plot(r * 1e10, probability_density,'LineWidth',2);
xlabel('Radius (Å)');
ylabel('Probability Density');
title(['Probability Density of ', atom,' Atom in ', orbital]);
grid on;

theta = linspace(0, pi, 50);
phi = linspace(0, 2*pi, 50);
[Theta,Phi] = meshgrid(theta, phi);

switch orbital
    case '1s'
        R = (Z^(3/2) / sqrt(pi)) * exp(-Z * r_min);
    case '2s'
        R = (Z^(3/2) / (4*sqrt(2))) * (2 - Z * r_min) .* exp(-Z * r_min / 2);
    case '2p'
        R = (Z^(3/2) / (4*sqrt(6))) * Z * r_min .* exp(-Z * r_min / 2);
    case '3s'
        R = (Z^(5/2) / (27*sqrt(6))) * (27 - Z * r_min) .* exp(-Z * r_min / 3);
    case '3p'
        R = (Z^(5/2) / (27*sqrt(6))) * Z * r_min .* exp(-Z * r_min / 3);
    case '3d'
        R = (Z^(7/2) / (81*sqrt(30))) * Z^2 * r_min.^2 .* exp(-Z * r_min / 4);
    case '4s'
        R = (Z^(7/2) / (64*sqrt(24))) * (64 - Z * r_min) .* exp(-Z * r_min /4);
end

X=R.*sin(Theta).*cos(Phi);
Y=R.*sin(Theta).*sin(Phi);
Z_coord=R.*cos(Theta);

figure;
surf(X,Y,Z_coord,'EdgeColor','none');
xlabel('X');
ylabel('Y');
zlabel('Z');
title(['3D Plot of ', atom,' Atom in ', orbital]);
axis equal;
view(30,30); 
colorbar;
