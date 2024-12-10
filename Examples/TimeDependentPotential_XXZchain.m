% Example: Dynamics of an XXZ Spin Chain in a Time-Varying Magnetic Field
%
% This script demonstrates the simulation of a 1D XXZ spin chain subject to 
% a parabolic magnetic confinement, whose strength decreases over time. 
% The evolution of the system is simulated using the iFluid library. 
%
% Key calculations performed:
% 1. Initialize the thermal state of the XXZ chain in the presence of an 
%    initial magnetic confinement.
% 2. Propagate the state using the Backwards Semi-Lagrangian (BSL) solver 
%    with Runge-Kutta 4th order time-stepping.
% 3. Evaluate the dynamics of local conserved charges (density and energy) 
%    and their associated currents over time.
% 4. Visualize the spatiotemporal evolution of charges and currents at 
%    selected time steps.
%
% This example illustrates the effect of a time-dependent external field on 
% an integrable system and demonstrates how to use iFluid to study dynamics 
% in the XXZ spin chain.


clear all; 

iFluid_root_path = 'C:\Users\FrederikMoller\Documents\iFluid';
add_iFluid_path(iFluid_root_path)


%% Define simulation parameters
% Set grid and simulation parameters
N           = 2^7;                              % Number of rapidity grid points             
M           = 2^7;                              % Number of spatial grid points
dt          = 0.01;                             % Time step length
Ntypes      = 3;                                % Number of quasi-particle types

% Physical boundaries for rapidity and space
rmax        = pi/2;                             % Maximum rapidity value
zmax        = 1.5;                              % Maximum spatial position
tmax        = 1;                                % Total simulation time

% Define rapidity and spatial grids
[r_grid, rw] = legzo(N, -rmax, rmax);           % Generate Gauss-Legendre grid for rapidity
z_grid      = linspace(-zmax, zmax, M);         % Uniform spatial grid
t_array     = linspace(0, tmax, tmax/dt + 1);   % Time array for simulation steps

%% Define physical couplings and temperature
% Magnetic field confinement as a time-varying parabolic potential
syms x t                                        % Symbolic variables for position and time
B           = -1 - (1 - tanh(3 * t)) * 10 * x.^2; % Parabolic field with decreasing strength
B_func      = matlabFunction(B);               % Convert to a MATLAB function
dBdt        = matlabFunction(diff(B, t));      % Time derivative of the magnetic field
dBdx        = matlabFunction(diff(B, x));      % Spatial derivative of the magnetic field

% Define couplings: B_func (field), dBdt (time derivative), and dBdx (spatial derivative)
couplings  = { B_func    , @(t, x) acosh(1.5)   ; % Magnetic field coupling   
               dBdt      , []                  ; % Time derivative of coupling  
               dBdx      , []                  }; % Spatial derivative of coupling

T       = 1;                                   % Temperature of the initial thermal state

%% Initialize state and solve dynamics
% Initialize the XXZ spin chain model
XXZ         = XXZchainModel(z_grid, r_grid, rw, couplings, Ntypes);

% Compute the initial thermal state for the given temperature
theta_init  = XXZ.calcThermalState(T);

% Define and initialize the solver for dynamics
Solver      = AdvectionSolver_BSL_RK4(XXZ, 'implicit', false); % RK4 solver
fill_t      = Solver.propagateFilling(theta_init, t_array);      % Propagate the state over time

% Update couplings for dynamics after time evolution
XXZ.setCouplings( {@(t, x) 0 , @(t, x) acosh(1.5)} );          % Set updated field

% Compute the conserved charges and currents over time
[n_t, j_t]  = XXZ.calcCharges([0 1 2], fill_t, t_array);      % Charges and currents for selected indices

%% ------------ Plot results -------------------
% Time steps to sample for visualization
t_sample = [1 31 61 101]; % Selected time steps

% Create a figure to plot the evolution of charges and currents
figure

% Plot charge density q_0
sax1 = subplot(2, 2, 1);
plot(z_grid, n_t(:, t_sample, 1), 'LineWidth', 1.5)
xlim([-zmax, zmax])
xlabel('x')
ylabel('q_0')

% Plot charge density q_2
sax2 = subplot(2, 2, 2);
plot(z_grid, n_t(:, t_sample, 3), 'LineWidth', 1.5)
xlim([-zmax, zmax])
xlabel('x')
ylabel('q_2')

% Plot current j_0
sax3 = subplot(2, 2, 3);
plot(z_grid, j_t(:, t_sample, 1), 'LineWidth', 1.5)
xlim([-zmax, zmax])
xlabel('x')
ylabel('j_0')

% Plot current j_2
sax4 = subplot(2, 2, 4);
plot(z_grid, j_t(:, t_sample, 3), 'LineWidth', 1.5)
xlim([-zmax, zmax])
xlabel('x')
ylabel('j_2')

% Add a legend and title
legend(strcat('t=', string(num2cell(t_array(t_sample)))))
sgtitle('Evolution of Charges and Currents')
