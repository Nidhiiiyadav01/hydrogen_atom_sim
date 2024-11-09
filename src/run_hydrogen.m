% Define parameters for hydrogen atom
Z = 1;                         % Atomic number of hydrogen
basisParams = [1.0];           % Gaussian exponents for hydrogen 1s orbital
orbitalIndices = [1, 0];       % Quantum numbers [n, l] for 1s orbital

% Create an atom instance for hydrogen
hydrogen = Atom(Z, basisParams, orbitalIndices);

% Build matrices and solve eigenproblem
hydrogen = hydrogen.build_matrices();
hydrogen = hydrogen.solve_eigenproblem();

% Display calculated energy eigenvalues
disp('Energy Eigenvalues for Hydrogen:');
disp(hydrogen.E);

% Visualize the 1s atomic orbital
grid_size = 50;     % Resolution of the grid
grid_range = 1.5;   % Range for x, y, z grid
orbital_index = 1;  % Index of the orbital to visualize

% Plot the atomic orbital density
hydrogen.visualize_orbital(orbital_index, grid_size, grid_range);
