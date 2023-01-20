clear all; close all;

% In this example dynamics of a partitioning protocol in the relativistic
% sinh-Gordon model are calculated.

% Add paths to iFluid directories
addpath(['..' filesep 'models' filesep])
addpath(['..' filesep 'solvers' filesep])
addpath(['..' filesep 'iFluid' filesep])
addpath(['..' filesep 'utils' filesep])

%% Define simulation parameters

N           = 2^7;                              % number of rapidity gridpoints             
M           = 2^7;                              % number of position gridpoints
dt          = 0.025;                            % time step length

rmax        = 3;                                % max rapidity
zmax        = 2;                                % max position
tmax        = 2;                                % max time

% use Gauss-Legendre quadrature for rapidity
[r_grid,rw] =legzo(N, -rmax, rmax);             % rapidity grid and quadrature weights
z_grid      = linspace(-zmax, zmax, M);
t_array     = linspace(0, tmax, tmax/dt+1);


%% Define physical couplings and temperature

alpha       = @(t,x) (1+ 0.5*tanh(2*t))/(8*pi + 2);
dadt        = @(t,x) sech(2*t).^2 /(8*pi + 2);

couplings   = { alpha  , @(t,x) 1 , @(t,x) 0 ;     % coupling    
                dadt   , []       , []       ;     % d/dt coupling  
                []     , []       , []       };    % d/dx coupling

T           = @(x) 1.5 + 0.25*tanh(50*x);          % position dependent temperature


%% Initialize solvers and model + initial state

shG         = sinhGordonModel(z_grid, r_grid, rw, couplings);

% enable extrapolation for homogeneous non-zero filling outside zmax
Solver      = AdvectionSolver_BSL_RK4(shG, 'implicit', false, 'extrapolate', true);

% balance density difference between two halves of system via offset
% chemical potential in initial state.
coup_init   = { alpha  , @(t,x) 1 , @(t,x) 0.289*tanh( - 50*x ) };

theta_init  = shG.calcThermalState(T, coup_init);


%% Solve dynamics and calculate expectation values
theta_t     = Solver.propagateTheta(theta_init, t_array);

[n_t, j_t]  = shG.calcCharges([0 1 2], theta_t, t_array); % (density, momentum and energy)

Psi_k_t     = shG.calcVertexExpval( 3, theta_t, t_array);


%% ------------ Plot results -------------------

sample_idx = [ 1    11    21    41    61    81 ];

figure
sax1 = subplot(2,1,1);
hold on
box on
for i = 1:length(sample_idx)
    plot(z_grid, Psi_k_t(:,3,sample_idx(i)),'LineWidth',1.5)
end
ylabel('\langle \Phi _{3} \rangle')
xlabel('position z')

sax2 = subplot(2,1,2);
hold on
box on
for i = 1:length(sample_idx)
    plot(z_grid, n_t(:,sample_idx(i),1),'LineWidth',1.5)
end
xlabel('position z')
ylabel('q_0')

legstr = [];
for i = 1:length(sample_idx)
    legstr{i} = ['t = ' num2str(t_array(sample_idx(i) )) ];
end

leg = legend(legstr,'NumColumns',2,'Location','Northeast');

legend('boxoff')
