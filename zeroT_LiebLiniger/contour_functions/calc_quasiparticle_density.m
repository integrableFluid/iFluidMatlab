function [rho, rapid_grids, rapid_weights] = calc_quasiparticle_density(fermi_points, c, N, Nmin)
    % =====================================================================
    % Purpose : Solve TBA equations to calculate quasiparticle distribution 
    %           for a given set of local Fermi rapidities.
    % Input :   fermi_points-- cell array, where each entry is a vector of
    %                           local Fermi rapidities
    %           c           -- Lieb-Liniger coupling constant
    %           N           -- number of points used to discretize rapidity
    %           Nmin        -- minimum number of points used to discretize
    %                           a Fermi sea.
    % Output:   rho         -- quasiparticle distribution for each set of 
    %                           Fermi points
    %           rapid_grids -- rapidity grid for each set of Fermi points
    %           rapid_weights -- quadrature weights for each Fermi points
    % =====================================================================

    if nargin < 4
        Nmin = ceil(N/10);
    end

    if ~iscell(fermi_points)
        % fermi-points should be cell
        fermi_points = {fermi_points};
    end


    dp          = @(rapid) ones(size(rapid));
    rho         = zeros(N, length(fermi_points));
    rapid_grids = zeros(N, length(fermi_points));
    rapid_weights = zeros(N, length(fermi_points));
    
    for i = 1:length(fermi_points)
        % discretize rapidity intervals between pairs for Fermi points
        [grid, weights] = create_adaptive_grids(fermi_points{i}(:), N, Nmin);

        % solve dressing equation
        dp_dr           = apply_dressing(dp, grid, weights, c);
        
        if isempty(dp_dr)
            dp_dr = zeros(N, 1);
            grid = zeros(N, 1);
            weights = zeros(N, 1);
        end
        
        % quasiparticle density is equal to density of states at T=0
        rho(:,i)        = dp_dr/(2*pi);   
        rapid_grids(:,i)= grid; 
        rapid_weights(:,i)= weights; 
    end
end