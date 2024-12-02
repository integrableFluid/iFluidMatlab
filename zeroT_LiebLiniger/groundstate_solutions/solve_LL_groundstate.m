function [g, gamma, e, Q] = solve_LL_groundstate(lambda, N)
    % Given the scalar lambda = c/K, where c is the coupling strength and
    % K is the Fermi momentum, solve the equations in the LL paper to
    % obtain the ground state functions.
    %
    % g is the quasi-particle distribution between -K and K
    % gamma is interaction strength
    % e is the scaled energy
    % Q is the Fermi momentum in units of linear density
    
    grid    = linspace(-1, 1, N)'; 
    delta   = grid(2) - grid(1);
    kernel  = -1/pi * lambda./(lambda^2 + (grid - grid').^2);
    source  = ones(size(grid))/(2*pi);

    % solve eq. (3.18) in original LL paper
    g       = (eye(length(grid)) + kernel.*delta)\source;

    % use eq. (3.20) to determine gamma
    gamma   = lambda/trapz(grid, g);

    % use eq. (3.19) to determine energy density
    e       = gamma^3/lambda^3 * trapz(grid, grid.^2 .* g);
    
    % use eq. (3.21) to calculate Fermi point K in units of density n
    Q       = gamma/lambda;
end