clear all; close all;


% Reproduces the results found in the paper:
% "Generalized hydrodynamics with space-time inhomogeneous interactions".
% (http://arxiv.org/abs/1906.01654)

addpath(['..' filesep 'models' filesep])
addpath(['..' filesep 'solvers' filesep])
addpath(['..' filesep 'iFluid' filesep])
addpath(['..' filesep 'utils' filesep])

%% Define simulation parameters

N           = 50; %2^6;
M           = 2^6;
dt          = 0.025;

kmax        = 3;
xmax        = 3;
tmax        = 12.5;

k_array     = linspace(-kmax, kmax, N);
kw = k_array(2) - k_array(1);
% [k_array,kw]=legzo(N, -kmax, kmax);
x_array     = linspace(-xmax, xmax, M);
t_array     = linspace(0, tmax, tmax/dt+1);


%% Define physical couplings and temperature
couplings1  = { @(t,x) 0.5 - 0.5*x.^2   , @(t,x) 0.3 + tanh(3*t)    ;        % couplings
                []                      , @(t,x) 3*sech(3*t)^2      ;        % d/dt
                @(t,x) -x               , []                        };       % d/dx
            
couplings2  = { @(t,x) 0.5 - 0.5*x.^4   , @(t,x) 0.3 + tanh(3*t)    ;        % couplings
                []                      , @(t,x) 3*sech(3*t)^2      ;        % d/dt
                @(t,x) -2*x.^3          , []                        };       % d/dx

            
T           = 0.5;


%% Initialize state and solve dynamics

% Harmonic potential
LLS1        = LiebLinigerModel(x_array, k_array, kw, couplings1);
TEST        = SecondOrderSolver(LLS1, []);
theta1_init = LLS1.calcThermalState(T);
theta1_t    = TEST.propagateTheta(theta1_init, t_array);
n1_t        = LLS1.calcCharges(theta1_t, 0, t_array);


% Anharmonic potential

LLS2        = LiebLinigerModel(x_array, k_array, kw, couplings2);
TEST2       = SecondOrderSolver(LLS2, []);
theta2_init = LLS2.calcThermalState(T);
theta2_t    = TEST2.propagateTheta(theta2_init, t_array);
n2_t        = LLS2.calcCharges(theta2_t, 0, t_array);



%% ------------ Plot results -------------------

% Initial state plot
figure
hold on
subplot(2,2,1)
imagesc(x_array,k_array, squeeze(double(theta1_init)) )
colormap(hot)
caxis([ 0 1])
set(gca,'YDir','normal') 
title('Harmonic')

subplot(2,2,2)
imagesc(x_array,k_array, squeeze(double(theta2_init)) )
colormap(hot)
caxis([ 0 1])
set(gca,'YDir','normal') 
title('Anharmonic')

subplot(2,2,3)
plot(x_array, n1_t(:,1))
xlabel('x')
ylabel('n')
ylim([0 1])
xlim([-xmax, xmax])

subplot(2,2,4)
plot(x_array, n2_t(:,1))
xlabel('x')
ylabel('n')
ylim([0 1])
xlim([-xmax, xmax])


% Stability plot
figure
hold on
box on
plot(t_array, sum(n1_t,1))
plot(t_array, sum(n2_t,1))
xlabel('t')
ylabel('N')
xlim([0 tmax])
legend('Harmonic','Anharmonic')

% find indices of specific times
plottimes       = [0, 0.1, 6, 12];
[~, t_idx]      = ismember(plottimes, t_array);



%% Harmonic Pontential
figure
for i = 1:4
    subplot(2,4, i)
    index = t_idx(i);
    imagesc(x_array,k_array, squeeze(double(theta1_t{index})) ) 
    colormap(hot)
    caxis([ 0 1])
    set(gca,'YDir','normal') 
    title(['t = ' num2str(t_array(index))])
end

subplot(2,4,[5 6])
hold on
box on
for i = 1:4
    plot(x_array, n1_t(:,t_idx(i)))
end
xlabel('x')
ylabel('n')
ylim([0 1])
xlim([-xmax, xmax])

subplot(2,4,[7 8])
box on
plot(t_array, n1_t(M/2, :))
xlabel('t')
ylabel('n')
ylim([0.5 1])
xlim([0 tmax])


%% Anharmonic Pontential
figure
for i = 1:4
    subplot(2,4, i)
    index = t_idx(i);
    imagesc(x_array,k_array, squeeze(double(theta2_t{index})) )
    colormap(hot)
    caxis([ 0 1])
    set(gca,'YDir','normal') 
    title(['t = ' num2str(t_array(index))])
end

subplot(2,4,[5 6])
hold on
box on
for i = 1:4
    plot(x_array, n2_t(:,t_idx(i)))
end
xlabel('x')
ylabel('n')
ylim([0 1])
xlim([-xmax, xmax])

subplot(2,4,[7 8])
box on
plot(t_array, n2_t(M/2, :))
xlabel('t')
ylabel('n')
ylim([0.5 1])
xlim([0 tmax])
xticks(0:2:12)