clear all; close all;

% In this example an XXZ chain is subject to a parabolic magnetic
% confinement, which decreases in strength over time.

% Add paths to iFluid directories
addpath(['..' filesep 'models' filesep])
addpath(['..' filesep 'solvers' filesep])
addpath(['..' filesep 'iFluid' filesep])
addpath(['..' filesep 'utils' filesep])

%% Define simulation parameters

N           = 2^7;                              % number of rapidity gridpoints             
M           = 2^7;                              % number of position gridpoints
dt          = 0.01;                             % time step length
Ntypes      = 3;                                % number of quasi-particle types

rmax        = pi/2;                             % max rapidity
zmax        = 1.5;                              % max position
tmax        = 1;                                % max time

% use Gauss-Legendre quadrature for rapidity
[r_grid,rw] =legzo(N, -rmax, rmax);             % rapidity grid and quadrature weights
z_grid      = linspace(-zmax, zmax, M);
t_array     = linspace(0, tmax, tmax/dt+1);


%% Define physical couplings and temperature
syms x t
B           = -1 - (1-tanh(3*t))*10*x.^2;
B_func      = matlabFunction( B );
dBdt        = matlabFunction( diff(B,t) );
dBdx        = matlabFunction( diff(B,x) );

couplings  = { B_func    , @(t,x) acosh(1.5)      ;     % coupling    
               dBdt      , []                    ;     % d/dt coupling  
               dBdx      , []                      };    % d/dx coupling
            

T       = 1;

%% Initialize state and solve dynamics

XXZ         = XXZchainModel(z_grid, r_grid, rw, couplings, Ntypes);
theta_init  = XXZ.calcThermalState(T);


Solver      = AdvectionSolver_BSL_RK4(XXZ, 'implicit', false);
theta_t     = Solver.propagateTheta(theta_init, t_array);
%%
XXZ.setCouplings( {@(t,x) 0 , @(t,x) acosh(1.5)} );

[n_t, j_t]  = XXZ.calcCharges([0 1 2], theta_t, t_array);


%% ------------ Plot results -------------------

t_sample = [1 31 61 101];

figure('Position',[500 200 600 350])

sax1 = subplot(2,2,1);
plot(z_grid,  n_t(:,t_sample,1),'Linewidth',1.5)
xlim([-zmax, zmax])
xlabel('x')
ylabel('q_0')

sax2 = subplot(2,2,2);
plot(z_grid,  n_t(:,t_sample,3),'Linewidth',1.5)
xlim([-zmax, zmax])
xlabel('x')
ylabel('q_2')

sax3 = subplot(2,2,3);
plot(z_grid,  j_t(:,t_sample,1),'Linewidth',1.5)
xlim([-zmax, zmax])
xlabel('x')
ylabel('j_0')

sax4 = subplot(2,2,4);
plot(z_grid,  j_t(:,t_sample,3),'Linewidth',1.5)
xlim([-zmax, zmax])
xlabel('x')
ylabel('j_2')

legend( strcat('t=',string(num2cell(t_array(t_sample)))) )
sgtitle('evolution of charges and currents')

