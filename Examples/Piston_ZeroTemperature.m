% Example: Zero-Temperature GHD and Dispersive Shock Waves in a Bose Gas
%
% This script demonstrates the zero-temperature formulation of Generalized 
% Hydrodynamics (GHD), where the state is described by Fermi seas parameterized 
% by a Fermi contour. The contour is discretized into Fermi points (x_i, K_i), 
% with filling equal to 1 inside the contour and 0 outside.
%
% Key features demonstrated in this example:
% 1. Ground state calculation: The Lieb-Liniger equation is solved to obtain 
%    the Fermi momentum of a homogeneous 1D Bose gas at zero temperature.
% 2. Dynamics of a piston:
%    - The piston accelerates the Fermi momenta, leading to the emergence of 
%      a second local Fermi sea, signaling the formation of a dispersive 
%      shock wave (DSW).
% 3. Comparison of predictions:
%    - **Whitham modulation theory**: Calculates the density profile within 
%      the shock wave, capturing its fine oscillatory structure.
%    - **Thermodynamic Bethe Ansatz (TBA)**: Coarser-scale predictions for 
%      the shock wave density, which do not capture rapid oscillations.
% 4. Tabulation of velocities: Highlights the advantage of pre-calculating 
%    and tabulating effective velocities for zero-temperature dynamics, where 
%    the range of Fermi momenta is often known in advance.
%
% This example provides a detailed exploration of zero-temperature GHD and 
% demonstrates the interplay between fine- and coarse-grained descriptions of 
% dispersive shock waves in integrable systems


clear all; 

iFluid_root_path = 'C:\Users\FrederikMoller\Documents\iFluid';
add_iFluid_path(iFluid_root_path)



%% Setup parameters 

% Physical system parameters
c           = 1;    % Lieb-Liniger coupling strength
dens0       = 300;  % background atomic density


% Piston parameters
vp_acc      = 4;    % piston acceleration
vp0         = 4;    % initial piston velocity
vp_func     = @(t) vp_acc*t + vp0;
xp_func     = @(t) 0.5*vp_acc*t^2 + vp0*t;


% Grid and simulation parameters
N           = 2^8;  % number of rapidity gridpoints for integral equations
M_init      = 15;   % number of Fermi points to discretize initial contour
tmax        = 2;    % duration of dynamics
Nsteps      = 50;   % number of time steps

t           = linspace(0, tmax, Nsteps);


%% Calculate Fermi momentum and correponding sound velocity

gamma0          = c/dens0;

% Load table of Lieb-Liniger ground state solutions; if this table does not
% exist, it can be generated using the script Generate_groundstate_lookuptable
LL_GS_tables    = load('LL_GS_tables.mat');

% Lookup Fermi momentum (in units of density) from table
K0              = interp1(LL_GS_tables.gamma, LL_GS_tables.K, gamma0, 'spline')*dens0;

% The sound velocity can be calculated in several different ways. Here we
% first discretize the Fermi sea rapidity interval [-K0, K0], then solve
% the dressing equation to compute the effective velocity, which we then
% evaluate at the Fermi edge K0.
[grid, weights] = create_adaptive_grids([-K0; K0], N, round(N/10) );
dp_dr           = apply_dressing(@(rapid) ones(size(rapid)), grid, weights, c);
de_dr           = apply_dressing(@(rapid) 2*rapid, grid, weights, c);
vsound          = de_dr(end)./dp_dr(end);


%% Calculate effective velocity table

% First, setup velocity tabulator
dK1_arr = linspace(2*K0, 2*K0+vp_func(tmax), 20); % array of 1st Fermi sea widths
dK2_arr = linspace(0, vp_func(tmax), 15); % array of 2nd Fermi sea widths
h_arr   = linspace(0, vp_func(tmax), 15); % array of gaps

% Initialize velocity tabulator and calculate lookup-tables for the
% effective velocity for states with 1 and 2 Fermi seas.
VTab    = effective_velocity_tabulator();

veff_map1 = VTab.generate_1FS_lookuptable(dK1_arr, c, N); 
veff_map2 = VTab.generate_2FS_lookuptable(dK1_arr, dK2_arr, h_arr, c, N);

VTab.set_1FS_lookuptable( veff_map1 );
VTab.set_2FS_lookuptable( veff_map2 );


% Create anonymous function, which calculates the effective velocity by
% looking it up in the table. Here the argument x refers to a Fermi contour
veff_func   = @(x) VTab.tabulate_veff_contour(x, c, N);


%% Setup initial Fermi contour and 

% Setup initial Fermi contour; this is neccesary when initial piston
% velocity is non-zero.
contour_init = setup_initial_piston_contour(vp_func(0), K0, vsound, M_init);

% Simulate piston dynamics; here we use a method that takes linear steps in
% time, however, there are also higher-order methods available.
tic
Gamma_t = propagate_contour_piston_linear(contour_init, t, K0, xp_func, vp_func, veff_func);
toc

%% Calculate density and velocity profile of dispersive shock wave using Whitham theory

idx = round(linspace(1,length(t), 6));
idx(1) = [];

figure
tiledlayout(2,length(idx),'TileSpacing', 'compact', 'padding','compact' )
ii = 1;

for i = idx
    % The piston solver only returns the Fermi points affected by the
    % piston; it is implicitly assumed that a background Fermi sea state
    % between [-K0, K0] exists. This function constructs the full Fermi
    % contour given the piston solution.
    contour_full = include_piston_background(Gamma_t{i}, K0, 150);

    % Automatically extract Riemann invariants and calculate DSW density
    [R, z_DSW]  = extract_Riemann_invariants({contour_full}, 1e3);
    [n_DSW, u_DSW] = calc_DSW_density(R, z_DSW, t(i), c);

    % Calculate density with TBA
    zmax        = 1.25*z_DSW(end);
    z_TBA       = linspace( xp_func(t(i)),  1.25*z_DSW(end), 100);
    [n_TBA, u_TBA] = calc_hydrodynamic_variables(contour_full, z_TBA, c, N, 20);


    % Plot Fermi contour
    nexttile(ii)
    box on
    plot_contours(contour_full, 'b', [255, 165, 0]/255)
    xlim([0, zmax])
    xlabel('x')
    ylabel('\theta')


    % plot densities
    nexttile(ii+length(idx))
    hold on
    box on
    plot(z_DSW, n_DSW/dens0)
    plot(z_TBA, n_TBA/dens0)
    xlim([0, zmax])
    xlabel('x')
    ylabel('n/n_0')

    ii = ii+1;
end




