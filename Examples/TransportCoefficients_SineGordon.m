% Example: Transport Coefficients in the Quantum Sine-Gordon Model
%
% This script calculates transport coefficients (Drude weights and Onsager
% matrix elements) for the quantum Sine-Gordon model with various coupling
% configurations and particle spectra. The computations include:
%
% 1. Generation of all possible particle spectrum configurations based on
%    a maximum number of particles and magnonic levels.
% 2. Setting up the Sine-Gordon model for each particle configuration and
%    solving the Thermodynamic Bethe Ansatz (TBA) at a specified temperature.
% 3. Calculating the Drude weight and Onsager matrix for topological charge
%    and momentum.
% 4. Plotting the transport coefficients as functions of the coupling
%    parameter $\beta^2 / 8\pi$ in logarithmic scale.
%
% The results showcase the transport properties for both reflectionless 
% and full Sine-Gordon models under different spectra.

clear all; 

% Add paths to iFluid directories
iFluid_root_path = ['..' filesep];
add_iFluid_path(iFluid_root_path)

% Enable LaTeX interpretation for plots
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultAxesFontSize', 12);

% Enable GPU computation for iFluid, if supported
global enable_iFluid_GPU
enable_iFluid_GPU = true;

%% Setup parameters

N_parts_max = 10; % Maximum number of particles allowed in the spectrum
max_levels  = 1;  % Maximum number of magnonic levels in the spectrum

% Generate all unique configurations of particles (reflectionless and general)
particles = generateUniqueArrays(N_parts_max-1, max_levels+1);

% Add pure breather configurations (reflectionless points)
particles = [particles, num2cell(1:(N_parts_max-2))];

% Define physical parameters
T           = 0.25;      % Temperature of the system
mu          = 0;         % Chemical potential
c_idx       = [0, 1];    % Conserved quantities (charge and momentum)

% Define rapidity grid for the TBA equations
N           = 501;                              % Number of rapidity grid points
rmax        = 10;                               % Maximum rapidity
r_grid      = linspace(-rmax, rmax, N)';        % Rapidities evenly spaced on grid

% Coupling configuration
couplings = {@(t, x) mu}; % Constant chemical potential across space and time

% Solver settings for the TBA equations
settings.tolerance = 1e-16;    % Convergence tolerance
settings.TBA_solver = 'Newton'; % Solver method for TBA
settings.maxcount = 25;        % Maximum iterations for TBA solver

%% Calculate transport coefficients for all particle configurations

% Initialize arrays for storing results
D_ab = zeros(length(c_idx), length(c_idx), length(particles)); % Drude weights
L_ab = zeros(length(c_idx), length(c_idx), length(particles)); % Onsager matrix
xi   = zeros(1, length(particles));                           % Coupling parameter

for i = 1:length(particles)
    % Extract particle spectrum for this configuration
    N_levels = length(particles{i});        % Number of levels in the spectrum
    particle_spectrum = particles{i};       % Current particle spectrum configuration

    % Initialize the Sine-Gordon model for this configuration
    if length(particle_spectrum) == 1
        % Reflectionless model (only breathers)
        SG = SineGordonModel_reflectionless(0, r_grid, r_grid(2) - r_grid(1), particle_spectrum, couplings, settings);
    else
        % Full Sine-Gordon model (solitons and breathers)
        SG = SineGordonModel(0, r_grid, r_grid(2) - r_grid(1), particle_spectrum, couplings, settings);
    end

    % Compute the coupling parameter xi = beta^2 / (8 pi)
    xi(i) = SG.getXi();

    % Solve the thermal TBA to find the filling functions
    fill = SG.calcThermalState(T);

    % Calculate transport coefficients
    D_ab(:, :, i) = SG.calcDrudeWeight(c_idx, fill); % Drude weights
    L_ab(:, :, i) = SG.calcOnsagerMatrix(c_idx, fill); % Onsager matrix
end

%% Plot the results

% Define symbols and colors for different spectra
symbols = {'o', 'd', '^', '*', 's', 'x', '<', '>', '+'};
colors  = lines(6); % Distinct colors for up to 6 levels

figure
tiledlayout(2, length(c_idx), 'TileSpacing', 'compact', 'Padding', 'compact');

% Plot Drude weights
for j = 1:length(c_idx)
    nexttile
    hold on; box on;
    for i = 1:length(particles)
        % Extract particle configuration and coupling parameter
        N_levels = length(particles{i});
        beta = xi(i) / (xi(i) + 1); % Dimensionless coupling parameter

        % Plot Drude weight for this configuration
        plot(beta, D_ab(j, j, i), 'LineStyle', 'none', ...
             'Marker', symbols{N_levels}, ...
             'Color', colors(N_levels, :), ...
             'MarkerSize', 3, ...
             'MarkerFaceColor', colors(N_levels, :));
    end
    xlabel('$\beta^2 /8 \pi$');
    ylabel('$\mathcal{D}$');
    set(gca, 'yscale', 'log'); % Logarithmic scale for transport coefficients
end

% Plot Onsager matrix elements
for j = 1:length(c_idx)
    nexttile
    hold on; box on;
    for i = 1:length(particles)
        % Extract particle configuration and coupling parameter
        N_levels = length(particles{i});
        beta = xi(i) / (xi(i) + 1); % Dimensionless coupling parameter

        % Plot Onsager matrix element for this configuration
        plot(beta, L_ab(j, j, i), 'LineStyle', 'none', ...
             'Marker', symbols{N_levels}, ...
             'Color', colors(N_levels, :), ...
             'MarkerSize', 3, ...
             'MarkerFaceColor', colors(N_levels, :));
    end
    xlabel('$\beta^2 /8 \pi$');
    ylabel('$\mathcal{L}$');
    set(gca, 'yscale', 'log'); % Logarithmic scale for transport coefficients
end

%% Helper function for generating unique particle arrays
function result = generateUniqueArrays(M, N)
    % Initialize result as an empty cell array
    result = {};

    % Recursive function to generate combinations
    function generateCombinations(currentComb, remainingSum, remainingSlots)
        if length(currentComb) >= 2 && remainingSum >= 0
            if length(currentComb) == 1 || (length(currentComb) > 1 && currentComb(end) >= 2)
                result{end + 1} = currentComb; % Add valid combination
            end
        end
        if remainingSlots == 0 || remainingSum < 0, return; end
        startVal = isempty(currentComb) * 0 + ~isempty(currentComb) * 1;
        for i = startVal:remainingSum
            generateCombinations([currentComb, i], remainingSum - i, remainingSlots - 1);
        end
    end

    % Start generating combinations
    generateCombinations([], M, N);
end
