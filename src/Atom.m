classdef Atom
    properties
        Z                % Atomic number of the atom
        basisParams      % Basis function parameters (exponents) for each orbital
        orbitalIndices   % Orbital indices and quantum numbers [n, l]
        H                % Hamiltonian matrix
        S                % Overlap matrix
        C                % Coefficients matrix
        E                % Energy eigenvalues
    end
    
    methods
        % Constructor to initialize the atom with atomic number and orbital basis parameters
        function obj = Atom(Z, basisParams, orbitalIndices)
            obj.Z = Z;
            obj.basisParams = basisParams;
            obj.orbitalIndices = orbitalIndices;
            obj.H = [];  % Initialize Hamiltonian matrix
            obj.S = [];  % Initialize Overlap matrix
            obj.C = [];  % Initialize coefficient matrix for orbitals
            obj.E = [];  % Initialize energy eigenvalues
        end
        
        % Function to calculate overlap integral for atomic orbitals
        function S = overlap_integral(~, alpha1, alpha2)
            S = (pi / (alpha1 + alpha2))^(3/2);
        end
        
        % Function to calculate Hamiltonian integral (kinetic + nuclear attraction)
        function H = hamiltonian_integral(obj, alpha1, alpha2)
            S = obj.overlap_integral(alpha1, alpha2);
            H = S * (-alpha1 * alpha2 / (alpha1 + alpha2) + obj.Z);
        end
        
        % Build Hamiltonian and Overlap matrices for atomic orbitals
        function obj = build_matrices(obj)
            numOrbitals = length(obj.basisParams);
            obj.H = zeros(numOrbitals, numOrbitals);
            obj.S = zeros(numOrbitals, numOrbitals);
            
            for i = 1:numOrbitals
                for j = 1:numOrbitals
                    obj.S(i, j) = obj.overlap_integral(obj.basisParams(i), obj.basisParams(j));
                    obj.H(i, j) = obj.hamiltonian_integral(obj.basisParams(i), obj.basisParams(j));
                end
            end
        end
        
        % Solve the generalized eigenvalue problem for H and S matrices
        function obj = solve_eigenproblem(obj)
            [obj.C, D] = eig(obj.H, obj.S);
            obj.E = diag(D); % Eigenvalues (orbital energies)
        end
        
        % Visualize atomic orbital based on quantum numbers (n, l)
        function visualize_orbital(obj, orbital_index, grid_size, grid_range)
            x = linspace(-grid_range, grid_range, grid_size);
            y = linspace(-grid_range, grid_range, grid_size);
            z = linspace(-grid_range, grid_size, grid_size);
            [X, Y, Z] = meshgrid(x, y, z);
            
            orbital_density = zeros(size(X));
            for i = 1:length(obj.basisParams)
                orbital_density = orbital_density + obj.C(i, orbital_index) * obj.gaussian_basis(X, Y, Z, obj.basisParams(i));
            end
            orbital_density = orbital_density .^ 2; % Probability density
            
            isovalue = max(orbital_density(:)) * 0.1;
            figure;
            h = patch(isosurface(X, Y, Z, orbital_density, isovalue));
            isonormals(X, Y, Z, orbital_density, h);
            set(h, 'FaceColor', 'cyan', 'EdgeColor', 'none');
            xlabel('x'); ylabel('y'); zlabel('z');
            orbital_name = obj.get_orbital_name(orbital_index);
            title(sprintf('Atomic Orbital %s Visualization', orbital_name));
            axis equal; view(3); camlight; lighting gouraud;
            colorbar;
        end
        
        % Get the orbital name based on its quantum numbers
        function orbital_name = get_orbital_name(obj, orbital_index)
            n = obj.orbitalIndices(orbital_index, 1);
            l = obj.orbitalIndices(orbital_index, 2);
            subshells = ['s', 'p', 'd', 'f', 'g', 'h'];
            if l + 1 <= length(subshells)
                orbital_name = sprintf('%d%s', n, subshells(l + 1));
            else
                orbital_name = sprintf('%d?', n); % Placeholder if l is out of bounds
            end
        end
    end
    
    methods (Static)
        % Gaussian basis function (s-type orbital)
        function val = gaussian_basis(x, y, z, alpha)
            r2 = x.^2 + y.^2 + z.^2;
            val = exp(-alpha * r2);
        end
    end
end
