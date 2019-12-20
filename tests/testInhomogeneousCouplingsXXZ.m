clear all; close all;


% Reproduces the results found in the paper:
% "Generalized hydrodynamics with space-time inhomogeneous interactions".
% (http://arxiv.org/abs/1906.01654)

addpath(['..' filesep 'all' filesep])

%% Define simulation parameters

N           = 2^6;
M           = 2^7;
Ntypes      = 2;
dt          = 0.0125;

kmax        = pi/2;
xmax        = 1;
tmax        = 0.3;

% k_array     = linspace(-kmax, kmax, N);
[k_array,kw]=legzo(N, -kmax, kmax);
x_array     = linspace(-xmax, xmax, M);
t_array     = linspace(0, tmax, tmax/dt+2);


T           = 0.25;

%% Define physical couplings
syms x t
coup_theta  = acosh(1.5 + 0.3*tanh(3*t).*sin( 4*pi*(x-t) ) );
Coup        = matlabFunction( coup_theta );
dCoup_dt    = matlabFunction( diff(coup_theta,t) );
dCoup_dx    = matlabFunction( diff(coup_theta,x) );


couplings   = {@(t,x) -1 - 8*x.^2   , Coup          ;       % couplings
               []                   , dCoup_dt      ;       % d/dt couplings
               @(t,x) -16*x         , dCoup_dx      };      % d/dx couplings


%% Initialize state and solve dynamics
XXZ         = XXZchainModel(x_array, k_array, kw, couplings, Ntypes);
theta_init  = XXZ.calcThermalState(T);

Solver2     = SecondOrderSolver(XXZ, []);
theta_t     = Solver2.propagateTheta(theta_init, t_array);


%% Calculate charges
sz_t        = XXZ.calcCharges(theta_t, 0, t_array);


couplings{1,1} = @(t,x) 0;
XXZ.setCouplings(couplings);
h_t         = XXZ.calcCharges(theta_t, 2, t_array) ; 


%% Plot initial state

figure
for i = 1:Ntypes
    subplot(1,Ntypes,i)
    imagesc(x_array,k_array,squeeze( theta_init(:,:,i,:,:) ))
    colormap(hot)
%     caxis([ 0 1])
    set(gca,'YDir','normal') 
end


figure
for i = 1:Ntypes
    subplot(1,Ntypes,i)
    surf(x_array,k_array,squeeze( theta_init(:,:,i,:,:) ))
%     colormap(hot)
%     caxis([ 0 1])
%     set(gca,'YDir','normal') 
end


figure
subplot(2,1,1)
plot(x_array, sz_t(:,1))
xlim([-0.7, 0.7])
ylim([0 0.5])
xlabel('x/L')
ylabel('Sz')

subplot(2,1,2)
plot(x_array, h_t(:,1))
xlim([-0.7, 0.7])
% ylim([-1 0])
xlabel('x/L')
ylabel('h')


%% ------------ Plot results -------------------

% find indices of specific times
plottimes       = [0, 0.2, 0.3];
[~, t_idx]      = ismembertol(plottimes, t_array, 1e-4);


figure
subplot(2,1,1)
plot(x_array, sz_t(:,t_idx))
xlim([-0.7, 0.7])
ylim([0 0.5])
xlabel('x/L')
ylabel('Sz')

subplot(2,1,2)
plot(x_array, h_t(:,t_idx))
xlim([-0.7, 0.7])
ylim([-1 0])
xlabel('x/L')
ylabel('h')


% %% Run movie
% figure
% for i = 1:length(t_array)
%     dist = theta_t{i};
%     
%     subplot(1,2,1)
%     imagesc(x_array,k_array,squeeze(dist(:,:,1,:,:)))
%     colormap(hot)
%     caxis([ 0 1])
%     ylabel('rapidity')
%     set(gca,'YDir','normal') 
%     
%     subplot(1,2,2)
%     imagesc(x_array,k_array,squeeze(dist(:,:,2,:,:)))
%     colormap(hot)
%     caxis([ 0 0.02])
%     ylabel('rapidity')
%     set(gca,'YDir','normal') 
%     
%     suptitle(['t = ' num2str(t_array(i))])
%     pause(0.05)
% end