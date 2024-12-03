function q = calc_charge_densities(charge_nr, fermi_points, c, N)
    % =====================================================================
    % Purpose : calulate densities of charges given local Fermi points.
    % Input :   charge_nr   -- index of charge (0 is number, 1 momentum, ...)
    %           fermi_points-- cell array, where each entry is a vector of
    %                           local Fermi rapidities
    %           c           -- Lieb-Liniger coupling constant
    %           N           -- number of points used to discretize rapidity
    % Output:   q           -- charge density for each set of Fermi points
    % =====================================================================

    if ~iscell(fermi_points)
        % make sure fermi_points is a cell array
        fermi_points = {fermi_points};
    end

    % bare charges of Lieb-Liniger model
    h   = @(rapid) rapid.^(charge_nr);   
    q   = zeros(length(fermi_points), 1);
    
    for i = 1:length(fermi_points)
        % discretize rapidity intervals between pairs for Fermi points
        [grid, weights] = create_adaptive_grids(fermi_points{i}(:), N, 5);
        % solve dressing equation
        h_dr            = apply_dressing(h, grid, weights, c);
        % calculate charge density
        q(i)            = sum(weights.*h_dr/(2*pi));        
    end
end