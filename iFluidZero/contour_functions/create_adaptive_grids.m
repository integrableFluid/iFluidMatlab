function [grid, weights, K_idx] = create_adaptive_grids(Fermi_points, N, N_min) 
    % =====================================================================
    % Purpose : discretize rapidity between pairs of Fermi points.
    % Input :   fermi_points-- vector of local Fermi_points
    %           N           -- number of points used to discretize rapidity
    %           Nmin        -- minimum number of points used to discretize
    %                           a Fermi sea.
    % Output:   grid        -- vector of rapidity gridpoitns
    %           weights     -- corresponding quadrature weights 
    %           K_idx       -- indices of grid of Fermi points
    % =====================================================================


    % Distribute gridpoints according to Fermi sea (FS) widths
    N_fs        = length(Fermi_points)/2; % number of FS
    endpoints   = reshape(Fermi_points, [2, N_fs]);
    widths      = endpoints(2,:) - endpoints(1,:); % widths of FS
    widths_frac = widths/sum(widths);
    widths_frac(isnan(widths_frac)) = 0;
    min_points  = N_min*ones(1,N_fs);
    min_points(widths == 0) = 1; % only single gridpoint for zero-width FS
    N_points    = max([min_points ; ceil(N*widths_frac)],[],1); % number of gridpoints per FS
    
    
    % Remove excess gridpoints (prioritise grids with many gridpoints)
    excess      = sum(N_points) - N; % number of excess gridpoints
    excess_fs   = round(excess*widths_frac);
    diff        = excess - sum(excess_fs);
    if diff ~= 0
        [~, idx] = max(excess_fs);
        excess_fs(idx) = excess_fs(idx) + diff;
    end
    
    N_points    = N_points - excess_fs;
    
    
    % Create grid
    grid        = [];
    weights     = [];
    K_idx       = zeros(2, N_fs);
    for i = 1:N_fs 
        p1          = endpoints(1,i);
        p2          = endpoints(2,i);
        
        if N_points(i) > 1
            [x,w]       = legzo(N_points(i), p1, p2); % Gauss-Legendre quadrature
        else
            x           = p1; % p1 and p2 should be the same in this case
            w           = 0; % zero integration weight
        end
            
        grid        = [grid; flip(x)'];
        weights     = [weights; flip(w)'];
        K_idx(:,i)  = [1 ; N_points(i)];
    end
    
    for i = 2:N_fs 
        K_idx(:,i)  = K_idx(:,i) + K_idx(2,i-1);
    end
end