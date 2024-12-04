function countour_full = include_piston_background(contour_piston, K0, L, M_edge)
    % =====================================================================
    % Purpose : Include background (Fermi sea state within [-K0, K0]) in
    %           countour solution of piston dynamics.
    % Input :   contour_piston -- Output contour of piston simulation,
    %                           containing only boosted Fermi points.
    %           K0      -- Background Fermi momentum/rapidity
    %           L       -- Max position of background // extend of system
    %           M_edge  -- (Optional) resolution of background
    % Output:   countour_full -- Full contour including background
    % =====================================================================

    if nargin < 4
        M_edge = 0;
    end

    % Detect whether piston component is fully disjoint from background
    is_2FS_state = contour_piston(end, 2) > K0;

    % first contour point should be at piston edge
    z_pist      = contour_piston(1,1); 
    K_pist      = contour_piston(1,2);

    if ~is_2FS_state % background is added as extra points in the contour
        % higher resolution back edge
        z_edge_b     = linspace(z_pist-eps, z_pist, M_edge+2)';
        K_edge_b     = linspace(-K0, K_pist, M_edge+2)';
    
        % higher resolution front edge
        z_edge_f     = linspace(L+eps, L, M_edge+2)';
        K_edge_f     = linspace(K0, -K0, M_edge+2)';
    
        countour_full= [ contour_piston;     
                        z_edge_f, K_edge_f
                        z_edge_b, K_edge_b;
                        ];

    else % background is added as a separate contour
        
        % setup top (piston) Fermi contour
        z_edge      = linspace(z_pist-eps, z_pist, M_edge+1)';
        K_edge      = linspace(contour_piston(end,2), K_pist, M_edge+1)';
        G1          = [ contour_piston; 
                        z_edge(1:end-1), K_edge(1:end-1);
                        ];
    
        % setup bottom (background) Fermi contour
        M_bg        = ceil(size(G1, 1)/4);
        x_bg        = linspace( z_pist, L, M_bg+2 )';
        K_bg        = linspace( -K0, K0, M_bg)';
        x_bg        = x_bg(2:end-1);
    
        G2          = [ x_bg, K0*ones(M_bg,1);
                        linspace(L+eps, L, M_bg)', flip(K_bg);
                        flip(x_bg), -K0*ones(M_bg,1);
                        linspace(z_pist-eps, z_pist, M_bg)', K_bg;
                        ];
    
        while size(G2,1) > size(G1,1)
            G2(end-1,:) = [];
        end

    
        countour_full= cat(3, G1, G2);    
    end

end