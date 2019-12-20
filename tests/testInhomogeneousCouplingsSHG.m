clear all; close all;


addpath(['..' filesep 'all' filesep])

%% Define simulation parameters

N           = 2^7;
M           = 2^7;
dt          = 0.025;

kmax        = 3;
xmax        = 1;
tmax        = 2;

% k_array     = linspace(-kmax, kmax, N);
% kw = k_array(2) - k_array(1);
[k_array,kw]=legzo(N, -kmax, kmax);
x_array     = linspace(-xmax, xmax, M);
t_array     = linspace(0, tmax, tmax/dt+1);

options.extrapFlag = true;

solopt.extrapFlag = true;

%% Define physical couplings and temperature
% couplings  = { @(t,x) 0 , @(t,x) 1 };
% T           = 1;

couplings  = { @(t,x) 2/(8*pi + 2)  , @(t,x) 1 - 0.5*tanh(5*t)      ;     % coupling    
               []                   , @(t,x) -2.5*sech(5*t).^2      ;     % d/dt coupling  
               []                   , []                            };    % d/dx coupling
            
% T       = @(x) 1;
T       = @(x) 1.5 + 0.25*tanh(2*pi*x);

%% Initialize state and solve dynamics

t_sample = floor(linspace( 1, length(t_array), 9  ));

% Harmonic potential
tic
shG         = sinhGordonModel(x_array, k_array, kw, couplings, options);
theta_init  = shG.calcThermalState(T);

Solver2     = SecondOrderSolver(shG,solopt);
theta_t     = Solver2.propagateTheta(theta_init, t_array);

n_t         = shG.calcCharges(theta_t, 0, t_array);
toc

%%
tic
Psi_k_t    = shG.calcVertexExpval( 5, theta_t(t_sample), t_array(t_sample));
toc


%% ------------ Plot results -------------------

% Initial state plot
figure
hold on
subplot(1,2,1)
imagesc(x_array,k_array, squeeze(double(theta_init)) )
colormap(hot)
caxis([ 0 1])
set(gca,'YDir','normal') 

subplot(1,2,2)
plot(x_array, n_t(:,1))
xlabel('x')
ylabel('n')
ylim([0 1])
xlim([-xmax, xmax])


% Stability plot
figure
hold on
box on
plot(t_array, sum(n_t,1))
xlabel('t')
ylabel('N')
xlim([0 tmax])


figure
for i = 1:length(Psi_k_t)
    subplot(3,3,i)
    plot(x_array, Psi_k_t{i})
end


%% Run movie
figure
for i = 1:length(t_array)
    theta = double(theta_t{i});
    
    subplot(2,1,1)
    imagesc(x_array,k_array,squeeze(theta))
    colormap(hot)
    caxis([ 0 0.5])
    ylabel('rapidity')
    set(gca,'YDir','normal') 
    
    subplot(2,1,2)
    plot(x_array, n_t(:,i) )
    xlim([ -xmax, xmax])
    ylim([0 max( n_t(:) )])
    xlabel('x')
    ylabel('density')
    
    suptitle(['t = ' num2str(t_array(i))])
    pause(0.05)
end
