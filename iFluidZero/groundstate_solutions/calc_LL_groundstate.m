function [rho, K] = calc_LL_groundstate(c, dens, N, LL_GS_table)
    % =====================================================================
    % Purpose : Calculate quasi-particle distribution of Lieb-Liniger model
    %           on rapidity grid of N points between -K and K.
    % Input :   c           -- Lieb-Liniger coupling constant
    %           dens        -- atomic density
    %           N           -- number of rapidity gridpoints
    %           LL_GS_table -- table of groundstate solutions (Fermi rapid
    %                           calculated for a number of interactions)
    % Output:   rho         -- quasiparticle distribution
    %           K           -- Fermi rapidity
    % =====================================================================    
    
    % calculate interaction strength
    gamma      = c./dens;
    
    % lookup Fermi momentum (in units of density) from table
    Q          = interp1(LL_GS_table.gamma, LL_GS_table.K, gamma, 'spline');
    
    % get Fermi points in proper units
    K          = Q.*dens; 

    % solve LL GS equations
    lambda     = c./K;
    
    rho        = zeros(N, length(lambda));
    for j = 1:length(lambda)
        rho(:,j)    = solve_LL_groundstate(lambda(j), N);
    end
end