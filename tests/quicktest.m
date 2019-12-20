clear all; close all;

set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultAxesFontSize',14);

addpath(['..' filesep 'models' filesep])
addpath(['..' filesep 'solvers' filesep])
addpath(['..' filesep 'iFluid' filesep])
addpath(['..' filesep 'utils' filesep])

%% Define simulation parameters

N           = 2^7;
M           = 2^7;
dt          = 0.025;

kmax        = 3;
xmax        = 3;
tmax        = 1;

k_array     = linspace(-kmax, kmax, N);
kw          = k_array(2) - k_array(1);
x_array     = linspace(-xmax, xmax, M);
t_array     = linspace(0, tmax, tmax/dt+1);





%% Define physical couplings and temperature
couplings   = { @(t,x) 0.5 - 0.5*x.^2   , @(t,x) 1    ;        % couplings
                []                      , @(t,x) 0    ;        % d/dt
                @(t,x) -x               , []                        };       % d/dx
            
            
T           = 0.5;


%% Initialize state and solve dynamics

coup_init = { @(t,x) 0.5 - 0.5*x.^2   , @(t,x) 1   };


LLS         = LiebLinigerModel(x_array, k_array, kw, couplings);
Solver2     = SecondOrderSolver(LLS, []);
theta_init  = LLS.calcThermalState(T, coup_init);
theta_t     = Solver2.propagateTheta(theta_init, t_array);

%% ====================== plot results ===========================

t_sample = [2 5 10 20 40];
g_n         = LLS.calcLocalCorrelator( 4, theta_t(t_sample), t_array(t_sample));

figure
plot(x_array, g_n)
