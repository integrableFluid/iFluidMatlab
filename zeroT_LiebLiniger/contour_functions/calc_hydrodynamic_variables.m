function [n, u] = calc_hydrodynamic_variables(Gamma, x, c, N, Nmin)
    % Calculate the hydrodynamic velocity and from that compute the phase

    fermi_points = find_contour_crossings(x, Gamma);
    [rho, rapid_grids, rapid_weights] = calc_quasiparticle_density(fermi_points, c, N, Nmin);
    
    % calculate density
    n   = sum(rho.*rapid_weights, 1);
    
    % calculate hydrodynamic velocity
    u   = sum(2*rho.*rapid_weights.*rapid_grids, 1)./n;
end