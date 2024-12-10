% Example: Traveling Wave State and Boundary Conditions in GHD
%
% This script sets up a traveling wave state by applying a boost to a 
% thermal equilibrium state. The dynamics of the system are then simulated 
% using the Backwards Semi-Lagrangian (BSL) solver.
%
% Key features demonstrated in this example:
% 1. Traveling wave initialization: A thermal state is boosted to create a 
%    traveling wave, providing an example of a non-equilibrium steady state.
% 2. Boundary conditions: Different types of boundary conditions are explored 
%    in the BSL solver to illustrate their impact on the simulation results.
% 3. Visualization: The filling function is plotted for different evolution
%    times for the various boundary conditions.



clear all; 

% Add paths to iFluid directories
iFluid_root_path = ['..' filesep];
add_iFluid_path(iFluid_root_path)

%% Define simulation parameters

% physical parameters
c           = 1;                                    % coupling strength
mu0         = 3;                                    % central chemical potential
T           = 1;                                    % temperature
L           = 4;                                   % length of system

% grid parameters
N           = 150;                                  % number of rapidity gridpoints              
M           = 200;                                  % number of position gridpoints
Nsteps      = 100;                                  % number of time steps

rmax        = 5;                                    % max rapidity
zmax        = L/2;                                  % max position
tmax        = 1;                                    % evolution duration

r_grid      = linspace(-rmax, rmax, N);             % rapidity grid
dr          = r_grid(2) - r_grid(1);                % quadrature weights
z_grid      = linspace(-zmax, zmax, M);             % position grid
t_array     = linspace(0, tmax, Nsteps+1);          % time steps


%% Setup Lieb-Liniger model with couplings

% Couplings are defined as a cell array containing anonymous functions with
% a time and space argument. Each column corresponds to separate couplings;
% for the Lieb-Liniger model, these are chemical potential and interactions
% strength. 
% For time- and space-independent couplings, there is no need to specify
% derivatives of couplings.
couplings   = { @(t,x) mu0   , @(t,x) c };
            
    
LL_model    = LiebLinigerModel(z_grid, r_grid, dr, couplings);


%% Calculate initial state

% Calculate thermal state by solving Yang-Yang equation for temperature T
fill_therm  = double( LL_model.calcThermalState(T) );

% Setup "boosted wave" state by locally boosting part of filling; for
% this one can use interpolation
boost_max   = 0.5;
boost_width = 0.2;
boost_pos   = -L/2 + 4*boost_width;
boost_profil= boost_max*exp(-(z_grid-boost_pos).^2 /(2*boost_width^2));

fill_boost  = zeros(size(fill_therm));
for i = 1:M
    fill_boost(:,i) = interp1(r_grid+boost_profil(i), fill_therm(:,i), r_grid, 'spline');
end

% Setup "travelling wave" state by removing boost at negative rapidities, 
% i.e. setting them equal to the thermal state
fill_wave   = fill_boost;
fill_wave(r_grid < 0,:) = fill_therm(r_grid < 0, :);

% Plot the different filling functions
figure
tiledlayout(3,1)

nexttile
imagesc(z_grid,r_grid,fill_therm)
set(gca,'YDir','normal')
colormap(bluewhitered)
colorbar
xlabel('x')
ylabel('\theta')
title('thermal state')

nexttile
imagesc(z_grid,r_grid,fill_boost)
set(gca,'YDir','normal')
colormap(bluewhitered)
colorbar
xlabel('x')
ylabel('\theta')
title('locally boosted state')

nexttile
imagesc(z_grid,r_grid,fill_wave)
set(gca,'YDir','normal')
colormap(bluewhitered)
colorbar
xlabel('x')
ylabel('\theta')
title('travelling wave state')


%% Solve Euler-scale dynamics for different boundary conditions

% For supported boundary conditions, see implementation of abstract
% AdvectionSolver_BSL class. The default value is 'none'.
boundary_conditions = {'open', 'none', 'periodic', 'reflective'};

% Setup explicit Runge-Kutta 4th order solver; the boundary conditions can
% be specified both at solver initialization or later.
solver_eul  = AdvectionSolver_BSL_RK4(LL_model, 'implicit', false);

fill_t_eul  = cell(length(boundary_conditions), length(t_array));

for i = 1:length(boundary_conditions)
    solver_eul.setBoundaryConditions( boundary_conditions{i} );

    % Simulate dynamics by solving GHD advection equation, using the 
    % bipartite thermal state as initial state. 
    tic
    fill_t_eul(i,:) = solver_eul.propagateFilling( fluidcell(fill_wave), t_array );
    toc
end


% plot evolution of filling function for different boundary conditions
N_plots     = 5; 
plot_idx    = round(linspace(1, Nsteps+1, N_plots));

figure
tiledlayout(length(boundary_conditions), N_plots, 'TileSpacing','tight','Padding','tight')

for j = 1:length(boundary_conditions)
    for i = plot_idx
        nexttile
        imagesc(z_grid, r_grid, double(fill_t_eul{j,i}))
        set(gca,'ydir', 'normal')
        xlabel('position')
        ylabel('rapidity')
        colormap(bluewhitered)
    end
    colorbar

end

sgtitle('Euler-scale dynamics with different boundary conditions')
