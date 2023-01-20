clear all; 

% In this example a thermal state is initialized in a cos-potential, which
% at time t=0 is quenched to a sin-potential.
% The evolution is calculated with periodic boundary conditions for both 
% Euler-scale and diffusive scale GHD.

% Add paths to iFluid directories
addpath(['..' filesep 'models' filesep])
addpath(['..' filesep 'solvers' filesep])
addpath(['..' filesep 'iFluid' filesep])
addpath(['..' filesep 'utils' filesep])

%% Define simulation parameters

% physical parameters
c           = 1;                                    % coupling strength
mu0         = 1;                                    % central chemical potential
T           = 3;                                    % temperature
L           = 4;                                    % length of system

V_init      = @(t,x) cos(2*pi*x/L);                 % initial potential
V_evol      = @(t,x) sin(2*pi*x/L);                 % evolution potential
dVdz_evol   = @(t,x) 2*pi/L*cos(2*pi*x/L);          % gradient of evolution potential


% grid parameters
N           = 100;                                  % number of rapidity gridpoints              
M           = 100;                                  % number of position gridpoints
Nsteps      = 300;                                  % number of time steps

rmax        = 5;                                    % max rapidity
zmax        = L/2;                                  % max position
tmax        = 2;                                    % evolution duration

r_grid      = linspace(-rmax, rmax, N);             % rapidity grid
dr          = r_grid(2) - r_grid(1);                % quadrature weights
z_grid      = linspace(-zmax, zmax, M);             % position grid
t_array     = linspace(0, tmax, Nsteps+1);          % time steps


%% Setup Lieb-Liniger model


couplings  = { @(t,x) mu0-V_evol(t,x)   , @(t,x) c ;    % couplings
               []                       , []       ;    % d/dt
               @(t,x) -dVdz_evol(t,x)   , []       };   % d/dx
            
    
LL_model   = LiebLinigerModel(z_grid, r_grid, dr, couplings);


%% Calculate filling function of initial thermal state

% setup initial potential
coup_init       = couplings;
coup_init{1,1}  = @(t,x) mu0 - V_init(t,x);

% calculate thermal state using the couplings specified in coup_init
fill_therm  = LL_model.calcThermalState(T, coup_init);

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



%% Solve Euler-scale dynamics 

% Setup explicit Runge-Kutta 4th order solver with periodic boundary
% conditions
solver_eul  = AdvectionSolver_BSL_RK4(LL_model, 'implicit', false, 'periodic_BC', true);

% Simulate dynamics by solving GHD advection equation, using the thermal 
% state as initial state. Note, in iFluid, the filling function is often
% referred to as "theta".
tic
fill_t_eul  = solver_eul.propagateTheta( fill_therm, t_array );
toc

% Calculate g2-function
g2_t_eul    = LL_model.calcLocalCorrelator(2, fill_t_eul, t_array);

%% Solve diffusive dynamics 

% Setup dissipation kernel, here diffusion, with the SETTLS scheme 
source      = Diffusion(LL_model, 'propagation_scheme', 'SETTLS');

% Setup solver with explicit Runge-Kutta 4th order scheme as base (periodic
% boundary conditions are automatically copied to dissipation kernel)
solver_diff = AdvectionSolver_BSL_RK4(LL_model, 'implicit', false, 'source', source, 'periodic_BC', true);

% Simulate dynamics by solving GHD disspative equation
tic
fill_t_diff = solver_diff.propagateTheta( fill_therm, t_array );
toc

% Calculate g2-function
g2_t_diff   = LL_model.calcLocalCorrelator(2, fill_t_diff, t_array);

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
end

for i = plot_idx
    nexttile
    imagesc(z_grid, r_grid, double(fill_t_diff{i}))
    set(gca,'ydir', 'normal')
    xlabel('position')
    ylabel('rapidity')
    colormap(bluewhitered)
end

sgtitle('evolution of filling function')

% plot density carperts
figure
tiledlayout(2,1)

nexttile
imagesc(t_array,z_grid, g2_t_eul)
xlabel('time')
ylabel('position')
colorbar

nexttile
imagesc(t_array,z_grid, g2_t_diff)
xlabel('time')
ylabel('position')
colorbar

sgtitle('evolution of g_2-function')
