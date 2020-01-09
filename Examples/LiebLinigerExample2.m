clear all; close all;

% In this example, the dynamics of a 1-dimensional Bose gas in an expanding
% box potential are simulated. As the box expands, the particles do "work" 
% on the sides of the box, thus lowering their kinetic energy leading to an
% effective cooling of the gas.


% Add paths to iFluid directories
addpath(['..' filesep 'models' filesep])
addpath(['..' filesep 'solvers' filesep])
addpath(['..' filesep 'iFluid' filesep])
addpath(['..' filesep 'utils' filesep])

%% Define simulation parameters

N           = 2^6;                              % number of rapidity gridpoints
M           = 2^7;                              % number of spatial gridpoints
dt          = 0.025;                            % length of timestep

rmax        = 5;                                % max rapidity
xmax        = 10;                               % max posistion
tmax        = 8;                                % max time

[r_array,rw]= legzo(N, -rmax, rmax);            % rapidity grid and weights
x_array     = linspace(-xmax, xmax, M);         % position grid
t_array     = linspace(0, tmax, tmax/dt+1);     % array of timesteps


%% Define physical couplings and temperature
syms x t                                        % symbolic chemical potential
mu          = 2 - 17*( -tanh(0.75*(x + 4 + 4*tanh(0.3*t) )) + ...
                        tanh(0.75*(x - 4 - 4*tanh(0.3*t))) + 2);

mu_func     = matlabFunction( mu );             % convert to anonymous function
dmu_dt      = matlabFunction( diff(mu,t) );     % take temporal derivative
dmu_dx      = matlabFunction( diff(mu,x) );     % take spatial derivative

couplings   = { mu_func     , @(t,x) 1    ;     % coupling    
                dmu_dt      , []          ;     % d/dt coupling  
                dmu_dx      , []          };    % d/dx coupling
            

T           = 4;                                % temperature


%% Initialize state and solve dynamics

% Initialize TBA for the Lieb-Liniger model
LLS         = LiebLinigerModel(x_array, r_array, rw, couplings);

% Intialize second order solver of the GHD equation
Solver2     = SecondOrderSolver(LLS, []); 

% Calculate thermal state at temperature T
theta_init  = LLS.calcThermalState(T);

% Propagate intial state according to the couplings
theta_t     = Solver2.propagateTheta(theta_init, t_array);


%% Calculate atomic density and kinetic energy density

% Set potential to zero, to only get kinetic contribution to energy
LLS.setCouplings( {@(t,x) 0 , @(t,x) 1} );

% Atomic density (charge index = 0), etomic density (charge index = 2)
q_t         = LLS.calcCharges([0 2], theta_t, t_array);


%%  ====================== plot results ===========================
figure
sub1 = subplot(1,2,1);
imagesc( x_array, r_array, theta_t{1}.getType(1, 'd') )
xlabel('x', 'FontSize', 14)
ylabel('rapidity', 'FontSize', 14)
caxis([0 1])
colormap(sub1,'hot')
annotation( gcf,'textbox',...
            sub1.Position,...
            'String',['t = ' num2str(t_array(1))],...
            'Color', 'white', ...
            'LineStyle','none',...
            'FontWeight','bold',...
            'FontSize',12,...
            'FitBoxToText','off');

sub2 = subplot(1,2,2);
imagesc( x_array, r_array, theta_t{end}.getType(1, 'd') )
xlabel('x', 'FontSize', 14)
ylabel('rapidity', 'FontSize', 14)
caxis([0 1])
colormap(sub2,'hot')
annotation( gcf,'textbox',...
            sub2.Position,...
            'String',['t = ' num2str(t_array(end))],...
            'Color', 'white', ...
            'LineStyle','none',...
            'FontWeight','bold',...
            'FontSize',12,...
            'FitBoxToText','off');


figure
for i = 1:6
    t_idx = ceil( (length(t_array)-1)/(6-1)*(i-1) + 1 );
   
    sax = subplot(6,1,i);
    hold on
    box on
    area(x_array, q_t(:, t_idx, 1));
    ylim([0 1.5])
    ylabel('n(x)', 'FontSize', 14)
    
    yyaxis right
    plot(x_array, 2 - mu_func(t_array(t_idx), x_array), 'LineWidth', 1.5 )
    ylim([0 40])
    ylabel('V(x)', 'FontSize', 14)
    
    if i < 6
        set(gca,'xticklabel',[])
    else
        xlabel('x', 'FontSize', 14)
    end
    
    annotation( gcf,'textbox',...
                sax.Position,...
                'String',['t = ' num2str(t_array(t_idx))],...
                'LineStyle','none',...
                'FontWeight','bold',...
                'FontSize',12,...
                'FitBoxToText','off');
end

figure
plot(t_array, trapz(x_array, q_t(:,:,2)))
xlabel('t', 'FontSize', 14)
ylabel('kinetic energy', 'FontSize', 14)
ylim([4 13])