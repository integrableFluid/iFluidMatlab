clear all; close all;

% In this example the dynamics of a quantum Newtons cradle are calculated.
% The cradle consist of two clouds of 1D Bose gases, which oscillate in a
% slightly anharmonic confinemt. 
% The evolution is calculated for both Euler-scale and diffusive scale GHD.

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

V_DW        = @(t,x) 0.15*x.^4 - 2*x.^2 + 0.1*x.^3; % double-well potential
dVdz_AHO    = @(t,x) 4*x.*exp(-2*x.^2/(12)^2);      % gradient of (an)harmonic oscillator potential


% grid parameters
N           = 150;                                  % number of rapidity gridpoints              
M           = 150;                                  % number of position gridpoints
Nsteps      = 500;                                  % number of time steps

rmax        = 8;                                    % max rapidity
zmax        = 6;                                    % max position
tmax        = 5;                                    % evolution duration

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


%% Solve diffusive dynamics 

% Setup dissipation kernel, here diffusion, with the SETTLS scheme 
source      = Diffusion(LL_model, 'propagation_scheme', 'SETTLS');

% Setup solver with explicit Runge-Kutta 4th order scheme as base 
solver_diff = AdvectionSolver_BSL_RK4(LL_model, 'implicit', false, 'source', source);

% Simulate dynamics by solving GHD disspative equation
tic
fill_t_diff = solver_diff.propagateTheta( fill_therm, t_array );
toc

% Calculate atomic density
n_t_diff    = LL_model.calcCharges(0, fill_t_diff, t_array);


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
imagesc(t_array,z_grid, n_t_eul)
xlabel('time')
ylabel('position')
colorbar

nexttile
imagesc(t_array,z_grid, n_t_diff)
xlabel('time')
ylabel('position')
colorbar

sgtitle('evolution of atomic density')