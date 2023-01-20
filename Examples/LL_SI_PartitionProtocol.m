

% Run partition protocol in external box trap (reflective boundary
% conditions). Each half of the system is initialized with the same atomic
% density, but different temperature.
% All quantities are in SI-units


clear all;



iFluid_path = ['..' filesep ];

addpath([iFluid_path filesep 'models' filesep])
addpath([iFluid_path filesep 'solvers' filesep])
addpath([iFluid_path filesep 'iFluid' filesep])
addpath([iFluid_path filesep 'utils' filesep])


%% Setup simulation parameters (these parameters can be tuned)

% physical parameters (in SI units)
L           = 80*1e-6;                          % boxlength
T_left      = 40*1e-9;                          % temperature of left half
T_right     = 80*1e-9;                          % temperature of right half
dens_init   = 75*1e6;                           % initial density
omegaT      = 2*pi*1.37*1e3;                    % transverse trapping frequency


% grid parameters
N           = 350;                              % number of rapidity gridpoints
M           = 100;                              % number of position gridpoints
Nsteps      = 200;                             % length of timestep

rmax        = 17.5*1e6;                         % max rapidity
zmax        = L/2;                              % max position 
tmax        = 50e-3;                           % max time



%% Setup grids and constant (these should not be tuned)

% grids 
r_grid      = linspace(-rmax, rmax, N);         % rapididty grid
z_grid      = linspace(-zmax, zmax, M);         % position grid
t_array     = linspace(0, tmax, Nsteps+1);      % array of timesteps
delta_r     = r_grid(2)-r_grid(1);


% constants
m_si        = 87*1.6605402e-27;                 % Rubidium-87 mass
hbar_si     = 1.054571726e-34;                  % reduced Planck constant
kB_si       = 1.38065e-23;                      % Boltzmann constant
as_si       = 5.2e-9;                           % Rubidium-87 s-wave scattering length

lperp_si    = sqrt(hbar_si/m_si/omegaT);        % transverse trapping width
c_si        = 2*as_si/lperp_si^2;               % coupling constant



%% Calculate initial state

% Here we use a trick to fit the density we want: Since each half of the
% system is supposed to be homogeneous, we only fit the atomic density for
% a single position "slice" and stack them next to each other to get the
% full phase-space distribution

couplings   = { @(t,x) 0*ones(1,M)      , @(t,x) c_si ;     % coupling    
                []                      , []          ;     % d/dt coupling  
                []                      , []          };    % d/dx coupling

% The "SI" LL-model is a wrapper class around the normal Lieb-Liniger
% model, which converts between SI-units and iFluid internal units
LLS         = LiebLinigerModel_SI(omegaT, 0, r_grid, delta_r, couplings); % set z_grid to 0

% calculate single-slice filling function for left and right side
[~,fill_left]  = LLS.fitDensity(T_left, dens_init, 0);
[~,fill_right] = LLS.fitDensity(T_right, dens_init, 0);


% stack the filling function slices next to each other to get full dist.
fill_init  = [ repmat(double(fill_left),1,M/2) , repmat(double(fill_right),1,M/2) ];
fill_init  = fluidcell(fill_init);


% plot the filling functions
figure
tiledlayout(2,1)

nexttile
hold on
box on
plot(r_grid, double(fill_left))
plot(r_grid, double(fill_right))
xlabel('rapidity')
ylabel('filling function')
legend('left','right')

nexttile
box on
imagesc(z_grid,r_grid, double(fill_init))
set(gca,'ydir','normal')
xlabel('position')
ylabel('rapidity')
colormap(bluewhitered)
colorbar


%% Simulate dynamics

% Re-initialize model, this time with correct z_grid
LLS         = LiebLinigerModel_SI(omegaT, z_grid, r_grid, delta_r, couplings); 

% Setup solver (implicit Runge-Kutta 2nd order with reflective BC)
solver      = AdvectionSolver_BSL_RK2(LLS, 'implicit', true, 'reflective_BC', true);

% run simulation (remember to convert time axis to internal units, as solver doesnt do this itself)
tic
theta_t     = solver.propagateTheta(fill_init, LLS.convert2TBA(t_array,'time') ); 
toc

% Calculate atomic density
n_t         = LLS.calcCharges(0, theta_t, t_array);

% Calculate energy density
e_t         = LLS.calcCharges(2, theta_t, t_array);

% Plot density carpets
figure
tiledlayout(2,1)

nexttile
imagesc(t_array,z_grid, n_t)
set(gca,'ydir','normal')
xlabel('time')
ylabel('position')
colorbar

nexttile
imagesc(t_array,z_grid, e_t)
set(gca,'ydir','normal')
xlabel('time')
ylabel('position')
colorbar

sgtitle('evolution of atomic density')
