% Example: Bipartition Protocol and Transition to Diffusive Dynamics
%
% This script demonstrates the calculation of particle and energy dynamics 
% in the bipartition protocol, where the system is initialized with an 
% offset chemical potential between the two halves. The evolution is studied 
% under both Euler-scale and diffusive-scale Generalized Hydrodynamics (GHD).
%
% Key calculations performed in this example:
% 1. Drude weight and Onsager matrix: Used to estimate the timescale at 
%    which diffusive dynamics become the dominant mode of transport for 
%    particles and energy.
% 2. GHD simulations: Euler-scale and diffusive-scale GHD are used to 
%    calculate particle and energy densities over time.
% 3. Visualization: Plots are generated to compare the scaling behavior 
%    of particle and energy densities, illustrating the transition from 
%    ballistic (Euler-scale) to diffusive dynamics.


clear all; 

% Add paths to iFluid directories
iFluid_root_path = ['..' filesep];
add_iFluid_path(iFluid_root_path)

%% Define simulation parameters

% physical parameters
c           = 1;                                    % coupling strength
mu0         = 2;                                    % central chemical potential
dmu         = 0.25;                                 % chemical potential offset
T           = 1;                                    % temperature
L           = 0.5;                                  % length of system

% grid parameters
N           = 150;                                  % number of rapidity gridpoints              
M           = 100;                                  % number of position gridpoints
Nsteps      = 100;                                  % number of time steps

rmax        = 5;                                    % max rapidity
zmax        = L/2;                                  % max position
tmax        = 0.025;                                % evolution duration

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
couplings  = { @(t,x) mu0   , @(t,x) c };
            
    
LL_model   = LiebLinigerModel(z_grid, r_grid, dr, couplings);


%% Calculate filling function of initial thermal state

% In bipartition, each half of the system is initialized with a different
% chemical potential.
% One can achieve this by having a Heaviside potential function, however, 
% this requires solving the Yang-Yang equation for many positions with
% identical chemical potential. 
% Alternatively, the thermal state at a single position in each half can be
% calculated, then duplicated to fill the entire grid.

% Left and right source terms of Yang-Yang eq.
source_L    = (LL_model.getBareEnergy(0, -L, r_grid', 1) - dmu)./T; 
source_R    = (LL_model.getBareEnergy(0, +L, r_grid', 1) + dmu)./T;      

% Solve Yang-Yang eq.
e_eff_L     = LL_model.calcEffectiveEnergy(source_L, 0, -L);
e_eff_R     = LL_model.calcEffectiveEnergy(source_R, 0, +L);

% Calculate filling
fill_L      = LL_model.calcFillingFraction(e_eff_L);
fill_R      = LL_model.calcFillingFraction(e_eff_R);

% Construct full filling function
fill_therm  = (z_grid < 0).*fill_L + (z_grid > 0).*fill_R; 

% Smoothen boundary slightly
fill_therm  = fluidcell( smoothdata(double(fill_therm), 2, 'gaussian', 5) );

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
ylabel('atomic density')

sgtitle('Initial state')



%% Calculate Drude weight and Onsager matrix to estimate diffusive time scale

% Calcuate thermal "background" state, i.e. thermal state for central
% chemical potential.
source0     = LL_model.getBareEnergy(0, 0, r_grid', 1)./T;      
e_eff0      = LL_model.calcEffectiveEnergy(source0, 0, 0);
fill0       = LL_model.calcFillingFraction(e_eff0);

% Calculate transport coefficients for density and energy (0th and 2nd charge)
c_idx       = [0, 2];
D_ij        = LL_model.calcDrudeWeight(c_idx, fill0)
L_ij        = LL_model.calcOnsagerMatrix(c_idx, fill0)

% The ratio of Onsager to Drude weight indicates the time-scale at which
% transport goes from diffusive to ballistic
t_ij        = L_ij./D_ij


%% Setup BSL solver and simulate Euler-scale and diffusive dynamics


% Setup explicit Runge-Kutta 4th order solver; for systems that extend
% beyond the specified position grid, use 'open' boundary conditions
Solver_Eul  = AdvectionSolver_BSL_RK4(LL_model, 'implicit', false, 'boundary_conditions', 'reflective');

tic
fill_t_Eul  = Solver_Eul.propagateFilling( fill_therm, t_array );
toc

% Setup dissipation kernel, here diffusion, with the SETTLS scheme 
source      = Diffusion(LL_model, 'propagation_scheme', 'SETTLS');

% Setup Adams-Moulton 4th order solver with source term
Solver_diff = AdvectionSolver_BSL_AM4(LL_model, 'boundary_conditions', 'reflective', 'source', source);

tic
fill_t_diff  = Solver_diff.propagateFilling( fill_therm, t_array );
toc


%% Calculate atomic and energy density

% Calculate atomic and energy density (0th and 2nd charges) for select
% evolution times.
t_eval  = t_array(1:10:end);
q_Eul   = LL_model.calcCharges([0, 2], fill_t_Eul(1:10:end), t_eval);
q_diff  = LL_model.calcCharges([0, 2], fill_t_diff(1:10:end), t_eval);

q_tot   = cat(4, q_Eul, q_diff);

%% Plot charges with ballistic and diffusive scaling
colors  = sky(length(t_eval));

fig_titles = {'atomic density', 'energy density'};
subfig_titles = {'Euler-scale', 'diffusive scale'};

for j = 1:2

    figure
    tile_handle = tiledlayout(3, 2,"TileSpacing","tight","Padding","tight");
    
    % plot profiles with no scaling
    for k = 1:2
        nexttile
        hold on
        box on
        for i = 1:length(t_eval)
            plot(z_grid, q_tot(:,i,j,k), 'Color', colors(i,:))
        end
        xlabel('x')
        ylabel('q')
        title(subfig_titles{k})  
    end

    % plot profiles with ballistic scaling
    for k = 1:2
        nexttile
        hold on
        box on
        for i = 1:length(t_eval)
            plot(z_grid/t_eval(i), q_tot(:,i,j,k), 'Color', colors(i,:))
        end
        xlabel('x/t')
        ylabel('q')
        xlim([-5, 5])
    end

    % plot profiles with diffusive scaling
    for k = 1:2
        nexttile
        hold on
        box on
        for i = 1:length(t_eval)
            plot(z_grid/sqrt(t_eval(i)), q_tot(:,i,j,k), 'Color', colors(i,:))
        end
        xlabel('x/\sqrt(t)')
        ylabel('q')
        xlim([-1, 1])
    end
    
    sgtitle(fig_titles{j})

    % Add a common colorbar
    colormap(colors);
    cb = colorbar();
    cb.Layout.Tile = 'east';
    cb.TickLabels = round(t_eval/t_ij(j,j), 3);
    cb.Label.String = 't/t_{crossover}'; % Label for the colorbar
end

