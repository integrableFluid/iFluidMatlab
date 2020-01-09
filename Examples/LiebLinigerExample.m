clear all; close all;

% In this example the dynamics of a quantum Newtons cradle are calculated.
% The cradle consist of two clouds of 1D Bose gases, which oscillate in a
% harmonic confinemt. As the clouds pass through each other, they interact.

% Add paths to iFluid directories
addpath(['..' filesep 'models' filesep])
addpath(['..' filesep 'solvers' filesep])
addpath(['..' filesep 'iFluid' filesep])
addpath(['..' filesep 'utils' filesep])

%% Define simulation parameters

N           = 2^7;
M           = 2^7;
dt          = 0.025;

rmax        = 13;
xmax        = 6;
tmax        = 8;

rap_array   = linspace(-rmax, rmax, N);
rap_w       = rap_array(2) - rap_array(1);
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


LLS        = LiebLinigerModel(x_array, rap_array, rap_w, couplings);
Solver2    = SecondOrderSolver(LLS, []);


offset      = 4;
coup_init   = { @(t,x) 2 - dw(t,x,offset)  , @(t,x) 1};
theta_init  = LLS.calcThermalState(T, coup_init);



[theta_t, u_t, w_t] = Solver2.propagateTheta(theta_init, t_array);
n_t        = LLS.calcCharges(0, theta_t, t_array);
rho_t      = LLS.transform2rho(theta_t, t_array);


%% ====================== plot density carpet ===========================

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
    imagesc(x_array, rap_array , rho_t{t_idx}.getType(1,'d') )
    set(gca,'YDir','normal')
    colormap(hot)
    caxis([0 0.3])
    
    
    if i == 1
        ylabel('\lambda')
    else
        set(gca,'yticklabel',[])
    end
    
end



%% ====================== plot characteristics ===========================

the_t   = theta_t{end}.getType(1,'d');
uj_t    = u_t{end}.getType(1,'d');
wj_t    = w_t{end}.getType(1,'d');
the_int = theta_init.getType(1,'d');

% interpolates theta(x) to u(t,x,lambda)
theta_u = interp2( x_array, rap_array, the_int, uj_t(:), wj_t(:), 'spline');
theta_u = reshape(theta_u, N, M);



figure
sax1 = subplot(2,2,1);
imagesc(rap_array, x_array, the_t )
set(gca,'YDir','normal') 
set(gca,'xticklabel',[])
caxis([0 1])

sax2 = subplot(2,2,2);
imagesc(rap_array, x_array,  theta_u )
set(gca,'YDir','normal')
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
caxis([0 1])


sax3 = subplot(2,2,3);
imagesc(rap_array, x_array, uj_t / xmax )
set(gca,'YDir','normal') 
caxis([-1.7 1.7])

sax4 = subplot(2,2,4);
imagesc(rap_array, x_array, wj_t / rmax)
set(gca,'YDir','normal') 
set(gca,'yticklabel',[])
caxis([-1.7 1.7])


colormap(sax1,'hot')
colormap(sax2,'hot')
colormap(sax3,'parula')
colormap(sax4,'parula')