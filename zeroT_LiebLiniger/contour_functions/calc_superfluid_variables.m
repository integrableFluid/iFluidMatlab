function [n, phi, u] = calc_superfluid_variables(fermi_points, x, c, N, Nmin)
    % =====================================================================
    % Purpose : calculate fluid quantities for given local Fermi points.
    % Input :   fermi_points-- cell array, where each entry is a vector of
    %                           local Fermi rapidities
    %           x           -- array of positions of Fermi points
    %           c           -- Lieb-Liniger coupling constant
    %           N           -- number of points used to discretize rapidity
    %           Nmin        -- minimum number of points used to discretize
    %                           a Fermi sea.
    % Output:   n           -- atomic density for each set of Fermi points
    %           phi         -- superfluid phase profile
    %           u           -- hydordynamic velocity profile
    % =====================================================================

    if nargin < 5
        Nmin = ceil(N/10);
    end
    
    % calculate quasi-particle density
    [rho, rapid_grids, rapid_weights] = calc_quasiparticle_density(fermi_points, c, N, Nmin);
    
    % calculate atomic density
    n   = sum(rho.*rapid_weights, 1);
    
    % calculate hydrodynamic velocity
    u   = sum(2*rho.*rapid_weights.*rapid_grids, 1)./n;

    % calculate superfluid phase (starting from 0)
    dx          = [0 reshape(diff(x), 1, [])];
    phi         = cumsum(dx.*u/2);
end