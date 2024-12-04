function countour_init = setup_initial_piston_contour(vp_init, K0, vsound, M)
    % =====================================================================
    % Purpose : Setup initial contour for piston problems.
    % Input :   vp_init -- Initial piston velocity
    %           K0      -- Background Fermi momentum/rapidity
    %           vsound  -- Sound velocity (effective velocity at K0)
    %           M       -- Number of points in initial contour
    % Output:   countour_init -- Initial contour
    % =====================================================================

    if vp_init > 0 % finite initial piston velocity
        % Create slanted edge in front of piston to avoid initially 
        % overlapping contour points. The rapidity of the edge depends on
        % whether the piston is above critial velocity.

        if vp_init < 2*vsound
            countour_init  = [flip(linspace(0, eps, M))', flip(linspace(K0, K0+vp_init, M))'];
        else
            countour_init  = [flip(linspace(0, eps, M))', flip(linspace(vp_init-K0, vp_init+K0, M))'];
        end
    
    else % piston starts at rest, i.e. it must be accelerating
        countour_init = [];
    end
    
end