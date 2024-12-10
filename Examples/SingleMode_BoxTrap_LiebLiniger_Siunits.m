% Example: Dynamics of the Lieb-Liniger Model in a Box Trap
%
% This script simulates the dynamics of a one-dimensional Bose gas described
% by the Lieb-Liniger model with reflective boundary conditions (box trap). 
% The system is initialized in a thermal state in a flat potential and then
% quenched to a spatially varying potential with a single cosine mode. 
% The script uses SI units, with seamless conversion to and from the internal
% units of the numerical solver via the SI extension of the Lieb-Liniger model.
%
% Key steps:
% 1. Define simulation parameters in SI units (temperature, chemical potential, etc.).
% 2. Initialize the system in a thermal state and fit the central chemical potential.
% 3. Set up a time-dependent quench potential and solve the dynamics using the GHD equations.
% 4. Compute and visualize the filling function and atomic density over time.

clear all; 

iFluid_root_path = 'C:\Users\FrederikMoller\Documents\iFluid';
add_iFluid_path(iFluid_root_path)


%% Setup simulation parameters
% Physical constants in SI units
m_si    = 87 * 1.6605402e-27;  % Mass of Rb-87 atom [kg]
hbar_si = 1.054571726e-34;     % Reduced Planck constant [J·s]
kB_si   = 1.38065e-23;         % Boltzmann constant [J/K]
as_si   = 5.2e-9;              % Scattering length [m]

% Numerical and system parameters
N        = 350;             % Number of rapidity grid points
Nsteps   = 200;             % Number of time steps
M        = 80;              % Number of spatial grid points
L        = 100e-6;          % Length of the box trap [m]
Natoms   = 5000;            % Total number of atoms in the system
rmax     = 15 * 1e6;        % Maximum rapidity [1/m]
tmax     = 100e-3;          % Total simulation time [s]
z_grid   = linspace(-L/2, L/2, M);  % Spatial grid [m]
r_grid   = linspace(-rmax, rmax, N); % Rapidity grid [1/m]
t_array  = linspace(0, tmax, Nsteps+1); % Time grid [s]

% Physical parameters
T        = 50e-9;                 % Temperature [K]
omega_perp = 2 * pi * 1.3e3;      % Transverse trapping frequency [Hz]
lperp    = sqrt(hbar_si / (m_si * omega_perp));  % Transverse confinement length [m]
c        = 2 * as_si / lperp^2;   % Interaction strength [1/m]

% Potential quench parameters
amp_mode    = 0.25;  % Amplitude of the quench mode as a fraction of chemical potential
nk_mode     = 2;     % Mode index (cosine spatial mode)
offset_mode = 0;     % Phase offset of the mode

%% Setup the Lieb-Liniger model and fit the initial chemical potential
% Initial potential (flat box trap)
V_init = @(t, x) 0;  

% Create the Lieb-Liniger model with initial couplings
couplings = { @(t,x) 0, @(t,x) c };  % Initial couplings (mu = 0, fixed c)
LL_model = LiebLinigerModel_SI(omega_perp, z_grid, r_grid, r_grid(2)-r_grid(1), couplings);

% Initial guess for the central chemical potential
mu_guess = 0.75 * hbar_si * omega_perp;

% Fit the chemical potential to match the desired atom number
mu0_fit = LL_model.fitAtomnumber(T, V_init, Natoms, mu_guess, true);

% Calculate the thermal state (filling function) at the fitted chemical potential
fill_init = LL_model.calcThermalState(T);

% Plot the initial filling function
figure;
box on;
plot(r_grid * 1e-6, fill_init(:, 1));  % Convert rapidity to [μm⁻¹] for the plot
xlabel('Rapidity $\theta$ [$\mu\mathrm{m}^{-1}$]', 'Interpreter', 'latex');
ylabel('Filling function $\vartheta$', 'Interpreter', 'latex');
title('Initial Filling Function', 'Interpreter', 'latex');

%% Solve dynamics after a quench
% Define time-dependent quench potential (single mode cosine potential)
dmudx = @(t, x) mu0_fit * amp_mode * 2 * pi * nk_mode / L * sin(offset_mode + 2 * pi * x * nk_mode / L);

% Update couplings for time-dependent dynamics
couplings = {
    @(t,x) mu0_fit,  @(t,x) c;          % Static couplings (mu and c)
    [],              [];                % Time derivatives of couplings (unused)
    @(t,x) dmudx(t,x), [];              % Spatial derivatives of couplings
};
LL_model.setCouplings(couplings);

% Set up the advection solver with reflective boundary conditions
Solver = AdvectionSolver_BSL_RK4(LL_model, ...
    'implicit', false, ...
    'boundary_conditions', 'reflective');

% Solve the GHD equations (convert time array to internal TBA units)
tic;
fill_t = Solver.propagateFilling(fill_init, LL_model.convert2TBA(t_array, 'time'));
toc;

% Compute atomic density at all times
n_t = LL_model.calcCharges(0, fill_t, t_array);

%% Plot the results
% Plot the evolution of the filling function (rapidity-space dynamics)
N_plots = 4;  % Number of snapshots to plot
plot_idx = round(linspace(1, Nsteps+1, N_plots));  % Indices of snapshots

figure;
tiledlayout(1, N_plots, 'TileSpacing', 'tight', 'Padding', 'tight');
for i = plot_idx
    nexttile;
    imagesc(z_grid * 1e6, r_grid * 1e-6, double(fill_t{i}));  % Convert to [μm]
    set(gca, 'YDir', 'normal');
    xlabel('$x$ [$\mu\mathrm{m}$]', 'Interpreter', 'latex');
    ylabel('$\theta$ [$\mu\mathrm{m}^{-1}$]', 'Interpreter', 'latex');
    colormap(bluewhitered);
    title(['$t = ' num2str(t_array(i) * 1e3) '\,\mathrm{ms}$'], 'Interpreter', 'latex');
end
sgtitle('Evolution of the Filling Function', 'Interpreter', 'latex');

% Plot the evolution of atomic density (real-space dynamics)
figure;
imagesc(t_array * 1e3, z_grid * 1e6, n_t);  % Convert to [ms] and [μm]
set(gca, 'YDir', 'normal');
xlabel('$t$ [$\mathrm{ms}$]', 'Interpreter', 'latex');
ylabel('$x$ [$\mu\mathrm{m}$]', 'Interpreter', 'latex');
colorbar;
title('Evolution of Atomic Density', 'Interpreter', 'latex');
