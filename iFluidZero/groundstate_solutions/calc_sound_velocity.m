function vsound = calc_sound_velocity(c, dens0, LL_GS_tables)
    % =====================================================================
    % Purpose : Calculate groundstate sound velocity.
    % Input :   c           -- Lieb-Liniger coupling constant
    %           dens        -- atomic density
    %           LL_GS_table -- table of groundstate solutions (Fermi rapid
    %                           calculated for a number of interactions)
    % Output:   vsound      -- sound velocity
    % =====================================================================    


    % calculate interaction strength
    gamma0      = c./dens0;
    
    % lookup Fermi momentum (in units of density) from table
    Q           = interp1(LL_GS_tables.gamma, LL_GS_tables.K, gamma0, 'spline');
    
    % get Fermi points in proper units
    K           = Q.*dens0; 
    
    % create grid between -K and K and solve dressing equation
    dp          = @(rapid) ones(size(rapid));
    de          = @(rapid) 2*rapid;
    [grid, w]   = create_adaptive_grids([-K ; K], 2^9, 5);
    dp_dr       = apply_dressing(dp, grid, w, c);
    de_dr       = apply_dressing(de, grid, w, c);
    
    % sounds velocity is given by effective velocity at Fermi point
    vsound      = de_dr(end)/dp_dr(end);
end