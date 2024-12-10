clear all; 

% Example: Dynamics in the Quantum Newton's Cradle Setup
%
% This script explores the dynamics of a system inspired by the quantum 
% Newton's cradle experiment. The initial state is prepared as a thermal 
% state in a tilted double-well potential, followed by a quench to a slightly 
% anharmonic potential.
%
% Key features demonstrated in this example:
% 1. Initial state preparation: Calculation of the quasiparticle distribution 
%    for a thermal state in a tilted double-well potential.
% 2. Dynamics simulation:
%    - Euler-scale GHD: Simulates the integrable dynamics.
%    - Quasi-1D GHD: Includes a collision integral to account for transverse 
%      excitations that break integrability.
% 3. Visualization: Comparison of the resulting dynamics for:
%    - Filling functions.
%    - Atomic density profiles.
%    - Populations of transverse states.
%
% This example provides insight into how anharmonicity and transverse 
% excitations affect the dynamics of a near-integrable system, bridging the 
% gap between idealized integrable models and real-world experiments.



clear all; 

% Add paths to iFluid directories
iFluid_root_path = ['..' filesep];
add_iFluid_path(iFluid_root_path)


%% Define simulation parameters

% physical parameters
c           = 1;                                    % coupling strength
mu0         = -3;                                    % central chemical potential
T           = 3;                                    % temperature

V_DW        = @(t,x) 0.15*x.^4 - 2*x.^2 + 0.05*x.^3; % double-well potential
dVdz_AHO    = @(t,x) 4*x.*exp(-2*x.^2/(12)^2);      % gradient of (an)harmonic oscillator potential


% grid parameters
N           = 100;                                  % number of rapidity gridpoints              
M           = 100;                                  % number of position gridpoints
Nsteps      = 500;                                  % number of time steps

rmax        = 8;                                    % max rapidity
zmax        = 6;                                    % max position
tmax        = 10;                                   % evolution duration

r_grid      = linspace(-rmax, rmax, N);             % rapidity grid
dr          = r_grid(2) - r_grid(1);                % quadrature weights
z_grid      = linspace(-zmax, zmax, M);             % position grid
t_array     = linspace(0, tmax, Nsteps+1);          % time steps


%% Setup Lieb-Liniger model

% note, the chemical potential and its derivative correspond to two
% different potentials, namely DW and HO. For calculating the thermal state
% (or energy density), the potential itself is used. FOr calculating the
% effective acceleration during dynamics, only its spatial derivative is
% used. 
% Since we wish to initialize the system in the DW and then evolve
% according to the HO, mixing the two potentials like this is OK.
couplings  = { @(t,x) mu0-V_DW(t,x)     , @(t,x) c ;    % couplings
               []                       , []       ;    % d/dt
               @(t,x) -dVdz_AHO(t,x)     , []       };   % d/dx
            
    
LL_model   = LiebLinigerModel(z_grid, r_grid, dr, couplings);


%% Calculate filling function of initial thermal state

% calculate thermal state using the couplings specified above and T
fill_therm  = LL_model.calcThermalState(T);

% calculate corresponding quasi-particle distribution
rho_therm   = LL_model.transform2rho(fill_therm);

% calculate corresponding atomic density (0th charge)
n_init      = LL_model.calcCharges(0, fill_therm, 0);

% plot the results. Note, the filling function is an object of the class
% "fluidcell". To plot it, cast it to a standard Matlab matrix via the
% function double().

figure
tiledlayout(3,1)

nexttile
imagesc(z_grid, r_grid, double(fill_therm))
set(gca,'ydir', 'normal')
xlabel('position')
ylabel('rapidity')
colormap(bluewhitered)
colorbar

nexttile
imagesc(z_grid, r_grid, double(rho_therm))
set(gca,'ydir', 'normal')
xlabel('position')
ylabel('rapidity')
colormap(bluewhitered)
colorbar

nexttile
plot(z_grid, n_init)
xlabel('position')
ylabel('density')


sgtitle('Initial thermal state in double-well potential')


%% Solve Euler-scale dynamics 

% Setup explicit Runge-Kutta 4th order solver
solver_eul  = AdvectionSolver_BSL_RK4(LL_model, 'implicit', false);

% Simulate dynamics by solving GHD advection equation, using the thermal 
% state as initial state. Note, in iFluid, the filling function is often
% referred to as "theta".
tic
fill_t_eul  = solver_eul.propagateTheta( fill_therm, t_array );
toc

% Calculate atomic density
n_t_eul     = LL_model.calcCharges(0, fill_t_eul, t_array);


%% Solve dynamics with quasi-1D collision integral

% Set length of transverse harmonic oscillator
lperp       = 0.35; 

% Setup dissipation kernel, here diffusion, with the SETTLS scheme 
source      = Quasi1D_CollisionIntegral(LL_model, lperp, 'propagation_scheme', 'endpoint');

% Setup solver with explicit Runge-Kutta 4th order scheme as base 
solver_q1D  = AdvectionSolver_BSL_RK4(LL_model, 'implicit', false, 'source', source);

% Simulate dynamics by solving GHD disspative equation
tic
[fill_t_q1D, trans_frac] = solver_q1D.propagateTheta( fill_therm, t_array );
toc

% Calculate atomic density
n_t_q1D     = LL_model.calcCharges(0, fill_t_q1D, t_array);


%% Plot the results 

% plot evolution of filling function
N_plots     = 4; 
plot_idx    = round(linspace(1, Nsteps+1, N_plots));

figure
tiledlayout(2, N_plots)

for i = plot_idx
    nexttile
    imagesc(z_grid, r_grid, double(fill_t_eul{i}))
    set(gca,'ydir', 'normal')
    xlabel('position')
    ylabel('rapidity')
    colormap(bluewhitered)
    yline(sqrt(2)/lperp, 'k--')
    yline(-sqrt(2)/lperp, 'k--')
end

for i = plot_idx
    nexttile
    imagesc(z_grid, r_grid, double(fill_t_q1D{i}))
    set(gca,'ydir', 'normal')
    xlabel('position')
    ylabel('rapidity')
    colormap(bluewhitered)
    yline(sqrt(2)/lperp, 'k--')
    yline(-sqrt(2)/lperp, 'k--')
end

sgtitle('Evolution of filling function')


% plot density carperts
figure
tiledlayout(2,1)

nexttile
imagesc(t_array,z_grid, n_t_eul)
xlabel('time')
ylabel('position')
colorbar

nexttile
imagesc(t_array,z_grid, n_t_q1D)
xlabel('time')
ylabel('position')
colorbar

sgtitle('Evolution of atomic density')


% plot fraction of population in transverse excited states
trans_frac = cell2mat(trans_frac(:));

figure
box on
hold on
plot(t_array, test(:,1))
plot(t_array, test(:,2))
xlabel('t')
ylabel('fractional transverse population')
