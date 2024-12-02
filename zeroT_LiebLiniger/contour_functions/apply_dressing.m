function Q_dr = apply_dressing(Q_func, rapid_grid, rapid_quad, c)
    % =====================================================================
    % Purpose : Solve Lieb-Liniger dressing equation for Fermi sea state.
    % Input :   Q_func      -- quantitity to be dressed, as anonnymous 
    %                           function of rapidity.
    %           rapid_grid  -- rapidity grid
    %           rapid_quad  -- rapidity integration quadrature weights
    %           c           -- Lieb-Liniger coupling constant
    % Output:   Q_dr        -- dressed quantity
    % =====================================================================
    Q           = Q_func(rapid_grid);
    kernel      = -1/pi * c./(c^2 + (rapid_grid - rapid_grid').^2);
    Q_dr        = (eye(length(rapid_grid)) + kernel.*rapid_quad')\Q;
end