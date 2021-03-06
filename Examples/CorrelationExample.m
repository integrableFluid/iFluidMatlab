clear all; close all;

% In this example, the spreading of Euler-scale two-point correlations is
% demonstrated. In this instance, density-density correlations in the
% Lieb-Liniger model are studied.

% Add paths to iFluid directories
addpath(['..' filesep 'models' filesep])
addpath(['..' filesep 'modules' filesep])
addpath(['..' filesep 'solvers' filesep])
addpath(['..' filesep 'iFluid' filesep])
addpath(['..' filesep 'utils' filesep])

%% Define simulation parameters

N           = 2^7;                              % number of rapidity gridpoints
M           = 2^8+1;                            % number of spatial gridpoints
dt          = 0.025;                            % length of timestep

rmax        = 7;                                % max rapidity
xmax        = 5;                                % max posistion
tmax        = 1;                                % max time

rap_array   = linspace(-rmax, rmax, N);         % rapidity grid
rap_w       = rap_array(2) - rap_array(1);      % rapidity quadrature weights
x_array     = linspace(-xmax, xmax, M);         % position grid
t_array     = linspace(0, tmax, tmax/dt+1);     % array of timesteps


%% Define physical couplings and temperature

% couplings are chemical potential and interaction strength
couplings  = { @(t,x) 0 , @(t,x) 1  ;           % couplings
               []       , []        ;           % d/dt couplings  
               []       , []        };          % d/dx couplings
                
T           = 0.5;                              % temperature                                 


%% Initialize state and solve dynamics

% Initialize TBA for the Lieb-Liniger model
TBA         = LiebLinigerModel(x_array, rap_array, rap_w, couplings);

% Intialize second order solver of the GHD equation
Solver2     = SecondOrderSolver(TBA, []);

% Set inital chemical potential and calculate thermal state
coup_init   = { @(t,x) 2 - 2*x.^2  , @(t,x) 1};
theta_init  = TBA.calcThermalState(T, coup_init);

% Propagate intial state according to the couplings
[theta_t, u_t] = Solver2.propagateTheta(theta_init, t_array);


%% Calculate correlations

% Initialize object for calculating correlations
iCorr       = CorrelationModule(TBA, []);

% Calculate form factors of density and associated current
calcFormFac = true;
[n,j,Vn,Vj] = TBA.calcCharges( 0, theta_t, t_array, calcFormFac );

% Calculate correlations at y = x_array((1+M)/2) and all times in t_array
[dir, indir]= iCorr.calc2PCorrelations( theta_t, u_t, t_array, (1+M)/2, Vn(2:end), Vj{1} );

%%  ====================== plot results ===========================

figure

subplot(2,1,1)
imagesc(t_array(2:end), x_array, squeeze(dir(:,:,1:end)))
cb = colorbar;
colormap(gca,blueredmap)
caxis([-1 1])
ylabel('Position x')
xlabel('Time t')

subplot(2,1,2)
imagesc(t_array(5:end), x_array, squeeze(indir(:,:,4:end)))
colorbar
colormap(gca,blueredmap)
ylabel('Position x')
xlabel('Time t')


%%

function map = blueredmap()

map  = [ 0         0    0.5000
         0    0.0076    0.5076
         0    0.0152    0.5152
         0    0.0227    0.5227
         0    0.0303    0.5303
         0    0.0379    0.5379
         0    0.0455    0.5455
         0    0.0530    0.5530
         0    0.0606    0.5606
         0    0.0682    0.5682
         0    0.0758    0.5758
         0    0.0833    0.5833
         0    0.0909    0.5909
         0    0.0985    0.5985
         0    0.1061    0.6061
         0    0.1136    0.6136
         0    0.1212    0.6212
         0    0.1288    0.6288
         0    0.1364    0.6364
         0    0.1439    0.6439
         0    0.1515    0.6515
         0    0.1591    0.6591
         0    0.1667    0.6667
         0    0.1742    0.6742
         0    0.1818    0.6818
         0    0.1894    0.6894
         0    0.1970    0.6970
         0    0.2045    0.7045
         0    0.2121    0.7121
         0    0.2197    0.7197
         0    0.2273    0.7273
         0    0.2348    0.7348
         0    0.2424    0.7424
         0    0.2500    0.7500
         0    0.2576    0.7576
         0    0.2652    0.7652
         0    0.2727    0.7727
         0    0.2803    0.7803
         0    0.2879    0.7879
         0    0.2955    0.7955
         0    0.3030    0.8030
         0    0.3106    0.8106
         0    0.3182    0.8182
         0    0.3258    0.8258
         0    0.3333    0.8333
         0    0.3409    0.8409
         0    0.3485    0.8485
         0    0.3561    0.8561
         0    0.3636    0.8636
         0    0.3712    0.8712
         0    0.3788    0.8788
         0    0.3864    0.8864
         0    0.3939    0.8939
         0    0.4015    0.9015
         0    0.4091    0.9091
         0    0.4167    0.9167
         0    0.4242    0.9242
         0    0.4318    0.9318
         0    0.4394    0.9394
         0    0.4470    0.9470
         0    0.4545    0.9545
         0    0.4621    0.9621
         0    0.4697    0.9697
         0    0.4773    0.9773
         0    0.4848    0.9848
         0    0.4924    0.9924
         0    0.5000    1.0000
    0.0152    0.5076    1.0000
    0.0303    0.5152    1.0000
    0.0455    0.5227    1.0000
    0.0606    0.5303    1.0000
    0.0758    0.5379    1.0000
    0.0909    0.5455    1.0000
    0.1061    0.5530    1.0000
    0.1212    0.5606    1.0000
    0.1364    0.5682    1.0000
    0.1515    0.5758    1.0000
    0.1667    0.5833    1.0000
    0.1818    0.5909    1.0000
    0.1970    0.5985    1.0000
    0.2121    0.6061    1.0000
    0.2273    0.6136    1.0000
    0.2424    0.6212    1.0000
    0.2576    0.6288    1.0000
    0.2727    0.6364    1.0000
    0.2879    0.6439    1.0000
    0.3030    0.6515    1.0000
    0.3182    0.6591    1.0000
    0.3333    0.6667    1.0000
    0.3485    0.6742    1.0000
    0.3636    0.6818    1.0000
    0.3788    0.6894    1.0000
    0.3939    0.6970    1.0000
    0.4091    0.7045    1.0000
    0.4242    0.7121    1.0000
    0.4394    0.7197    1.0000
    0.4545    0.7273    1.0000
    0.4697    0.7348    1.0000
    0.4848    0.7424    1.0000
    0.5000    0.7500    1.0000
    0.5152    0.7576    1.0000
    0.5303    0.7652    1.0000
    0.5455    0.7727    1.0000
    0.5606    0.7803    1.0000
    0.5758    0.7879    1.0000
    0.5909    0.7955    1.0000
    0.6061    0.8030    1.0000
    0.6212    0.8106    1.0000
    0.6364    0.8182    1.0000
    0.6515    0.8258    1.0000
    0.6667    0.8333    1.0000
    0.6818    0.8409    1.0000
    0.6970    0.8485    1.0000
    0.7121    0.8561    1.0000
    0.7273    0.8636    1.0000
    0.7424    0.8712    1.0000
    0.7576    0.8788    1.0000
    0.7727    0.8864    1.0000
    0.7879    0.8939    1.0000
    0.8030    0.9015    1.0000
    0.8182    0.9091    1.0000
    0.8333    0.9167    1.0000
    0.8485    0.9242    1.0000
    0.8636    0.9318    1.0000
    0.8788    0.9394    1.0000
    0.8939    0.9470    1.0000
    0.9091    0.9545    1.0000
    0.9242    0.9621    1.0000
    0.9394    0.9697    1.0000
    0.9545    0.9773    1.0000
    0.9697    0.9848    1.0000
    0.9848    0.9924    1.0000
    1.0000    1.0000    1.0000
    1.0000    1.0000    1.0000
    1.0000    0.9836    0.9836
    1.0000    0.9672    0.9672
    1.0000    0.9508    0.9508
    1.0000    0.9344    0.9344
    1.0000    0.9180    0.9180
    1.0000    0.9016    0.9016
    1.0000    0.8852    0.8852
    1.0000    0.8689    0.8689
    1.0000    0.8525    0.8525
    1.0000    0.8361    0.8361
    1.0000    0.8197    0.8197
    1.0000    0.8033    0.8033
    1.0000    0.7869    0.7869
    1.0000    0.7705    0.7705
    1.0000    0.7541    0.7541
    1.0000    0.7377    0.7377
    1.0000    0.7213    0.7213
    1.0000    0.7049    0.7049
    1.0000    0.6885    0.6885
    1.0000    0.6721    0.6721
    1.0000    0.6557    0.6557
    1.0000    0.6393    0.6393
    1.0000    0.6230    0.6230
    1.0000    0.6066    0.6066
    1.0000    0.5902    0.5902
    1.0000    0.5738    0.5738
    1.0000    0.5574    0.5574
    1.0000    0.5410    0.5410
    1.0000    0.5246    0.5246
    1.0000    0.5082    0.5082
    1.0000    0.4918    0.4918
    1.0000    0.4754    0.4754
    1.0000    0.4590    0.4590
    1.0000    0.4426    0.4426
    1.0000    0.4262    0.4262
    1.0000    0.4098    0.4098
    1.0000    0.3934    0.3934
    1.0000    0.3770    0.3770
    1.0000    0.3607    0.3607
    1.0000    0.3443    0.3443
    1.0000    0.3279    0.3279
    1.0000    0.3115    0.3115
    1.0000    0.2951    0.2951
    1.0000    0.2787    0.2787
    1.0000    0.2623    0.2623
    1.0000    0.2459    0.2459
    1.0000    0.2295    0.2295
    1.0000    0.2131    0.2131
    1.0000    0.1967    0.1967
    1.0000    0.1803    0.1803
    1.0000    0.1639    0.1639
    1.0000    0.1475    0.1475
    1.0000    0.1311    0.1311
    1.0000    0.1148    0.1148
    1.0000    0.0984    0.0984
    1.0000    0.0820    0.0820
    1.0000    0.0656    0.0656
    1.0000    0.0492    0.0492
    1.0000    0.0328    0.0328
    1.0000    0.0164    0.0164
    1.0000         0         0
    0.9918         0         0
    0.9836         0         0
    0.9754         0         0
    0.9672         0         0
    0.9590         0         0
    0.9508         0         0
    0.9426         0         0
    0.9344         0         0
    0.9262         0         0
    0.9180         0         0
    0.9098         0         0
    0.9016         0         0
    0.8934         0         0
    0.8852         0         0
    0.8770         0         0
    0.8689         0         0
    0.8607         0         0
    0.8525         0         0
    0.8443         0         0
    0.8361         0         0
    0.8279         0         0
    0.8197         0         0
    0.8115         0         0
    0.8033         0         0
    0.7951         0         0
    0.7869         0         0
    0.7787         0         0
    0.7705         0         0
    0.7623         0         0
    0.7541         0         0
    0.7459         0         0
    0.7377         0         0
    0.7295         0         0
    0.7213         0         0
    0.7131         0         0
    0.7049         0         0
    0.6967         0         0
    0.6885         0         0
    0.6803         0         0
    0.6721         0         0
    0.6639         0         0
    0.6557         0         0
    0.6475         0         0
    0.6393         0         0
    0.6311         0         0
    0.6230         0         0
    0.6148         0         0
    0.6066         0         0
    0.5984         0         0
    0.5902         0         0
    0.5820         0         0
    0.5738         0         0
    0.5656         0         0
    0.5574         0         0
    0.5492         0         0
    0.5410         0         0
    0.5328         0         0
    0.5246         0         0
    0.5164         0         0
    0.5082         0         0
    0.5000         0         0 ];
end