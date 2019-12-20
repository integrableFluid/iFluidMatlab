clear all; close all;


% Add paths to iFluid directories
addpath(['..' filesep 'models' filesep])
addpath(['..' filesep 'solvers' filesep])
addpath(['..' filesep 'iFluid' filesep])
addpath(['..' filesep 'utils' filesep])

%% Define simulation parameters

N           = 2^7 ;
M           = 2^7;
dt          = 0.01;
Ntypes      = 3;

kmax        = pi/2;
xmax        = 1.5;
tmax        = 1;

[k_array,kw]=legzo(N, -kmax, kmax);
x_array     = linspace(-xmax, xmax, M);
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

XXZ         = XXZchainModel(x_array, k_array, kw, couplings, Ntypes);
theta_init  = XXZ.calcThermalState(T);


Solver2     = SecondOrderSolver(XXZ, []);
theta_t     = Solver2.propagateTheta(theta_init, t_array);
%%
XXZ.setCouplings( {@(t,x) 0 , @(t,x) acosh(1.5)} );

[n_t, j_t]  = XXZ.calcCharges([0 1 2], theta_t, t_array);


%% ------------ Plot results -------------------

t_sample = [1 31 61 101];

figure('Position',[500 200 600 350])

sax1 = subplot(2,2,1);
plot(x_array,  n_t(:,t_sample,1),'Linewidth',1.5)
xlim([-xmax, xmax])
xlabel('x')
ylabel('q_0')

sax2 = subplot(2,2,2);
plot(x_array,  n_t(:,t_sample,3),'Linewidth',1.5)
xlim([-xmax, xmax])
xlabel('x')
ylabel('q_2')

sax3 = subplot(2,2,3);
plot(x_array,  j_t(:,t_sample,1),'Linewidth',1.5)
xlim([-xmax, xmax])
xlabel('x')
ylabel('j_0')

sax4 = subplot(2,2,4);
plot(x_array,  j_t(:,t_sample,3),'Linewidth',1.5)
xlim([-xmax, xmax])
xlabel('x')
ylabel('j_2')


