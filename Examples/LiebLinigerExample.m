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

kmax        = 13;
xmax        = 6;
tmax        = 8;

k_array     = linspace(-kmax, kmax, N);
kw          = k_array(2) - k_array(1);
x_array     = linspace(-xmax, xmax, M);
t_array     = linspace(0, tmax, tmax/dt+1);


%% Define physical couplings and temperature
sw    = @(t,x) 4*x.^2;
dw    = @(t,x,a) heaviside( -(x - a/2) ).*sw(t,x) + heaviside( (x - a/2) ).*sw(t,x-a);


couplings  = { @(t,x) 2 - 4*x.^2    , @(t,x) 1      ;        % couplings
               []                   , []            ;        % d/dt
               @(t,x) -8.*x         , []            };       % d/dx
            
    
T           = 3;


%% Initialize state and solve dynamics


LLS        = LiebLinigerModel(x_array, k_array, kw, couplings);
Solver2    = SecondOrderSolver(LLS, []);


offset      = 3;
coup_init   = { @(t,x) 2 - dw(t,x,offset)  , @(t,x) 1};
theta_init  = LLS.calcThermalState(T, coup_init);



[theta_t, u_t, w_t] = Solver2.propagateTheta(theta_init, t_array);
n_t        = LLS.calcCharges(0, theta_t, t_array);
rho_t      = LLS.transform2rho(theta_t, t_array);


%% ====================== plot results ===========================

figure

carp = subplot(2,6,1:6);
imagesc(t_array, x_array, n_t)
xlabel('Evolution time t')
ylabel('x')
caxis([0 2])

xticks([0 1.6 3.2 4.8 6.4 8])

for i = 1:6
    t_idx = ceil( (length(t_array)-1)/(6-1)*(i-1) + 1 );

    
    % plot quasiparticle distribution
    sax = subplot(2,6,i+6);
    imagesc(x_array, k_array , rho_t{t_idx}.getType(1,'d') )
    set(gca,'YDir','normal')
    colormap(hot)
    caxis([0 0.3])
    
    
    if i == 1
        ylabel('\lambda')
    else
        set(gca,'yticklabel',[])
    end
    
end


