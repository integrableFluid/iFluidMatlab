function [Gamma_t, t] = propagate_contour_ode45(Gamma, tspan, veff_func, ode_options)
    % =====================================================================
    % Purpose : Propagate Fermi contour with Matlab ODE45 solver.
    % Input :   Gamma       -- initial contour with each row = (z, rapid)
    %           tspan       -- ODE tspan (see Matlab documentation)
    %           veff_func   -- anonymous function to calculate effective
    %                           velocity (and acceleration) for contour
    %           ode_options -- ODE options (see Matlab documentation)
    % Output:   Gamma_t     -- cell array with contour for each time
    %           t           -- array of correposnding times
    % =====================================================================

    Nc      = size(Gamma, 3);
    Np      = size(Gamma, 1);
    
    y_init  = Gamma(:);

    if nargin < 4
        [t, y]  = ode45(@fp_derivs, tspan, y_init);
    else
        [t, y]  = ode45(@fp_derivs, tspan, y_init, ode_options);
    end
    
    Gamma_t = cell(1, length(t));
    for i = 1:length(t)
        Gamma_t{i} = reshape(y(i,:)', [Np, 2, Nc]);
    end
        
    function dydt = fp_derivs(t, y)
        G           = reshape(y, [Np, 2, Nc]);
        veff        = veff_func(G, 1:size(G,1));            
        dydt        = veff(:);
    end
end