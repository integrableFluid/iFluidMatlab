function contour = generate_contours(density, boost, c, z_grid, LL_GS_tables)
    % =====================================================================
    % Purpose : Generate Fermi contours for given set of density and boost.
    % Input :   density     -- density profiles (as rows)
    %           boost       -- boost profiles (as rows)
    %           c           -- Lieb-Liniger coupling constant
    %           z_grid      -- array of positions correponding to profiles
    %           LL_GS_tables-- table of groundstate LiebLiniger solutions
    % Output:   contour     -- 3D-array of (clockwise ordered) contours
    %                           each row is a set of points (z,rapid).
    % =====================================================================
    
    N_samples   = size(density, 1);
    M           = size(density, 2);
    
    % calculate interaction strength
    gamma       = c./reshape(density, 1, []);
    
    % lookup Fermi momentum (in units of density) from table
    Q           = interp1(LL_GS_tables.gamma, LL_GS_tables.K, gamma, 'spline');
    
    % reshape and get Fermi points in proper units
    Q           = reshape(Q, N_samples, M);
    K           = Q.*density; 
    
    % at each point in space z, the contour has two Fermi points shifted by
    % the boost, i.e. -K+\zeta and +K+\zeta    
    K_minus     = -K + boost;
    K_plus      = K + boost;
    
    % permute indices to have space 1st index and samples to 3rd index.
    % second index will be (z,\theta)-pairs
    K_minus     = permute(K_minus, [2 3 1]);
    K_plus      = permute(K_plus, [2 3 1]);
    
    % setup contours
    contour(:,1,:)  = repmat( [z_grid(:) ; flipud(z_grid(:))], 1,1,N_samples );
    contour(:,2,:)  = [ K_plus; flipud(K_minus) ];
end