function cross_list = find_contour_crossings(x, Gamma)
    % Find all intersections between vertical lines at points x and the set
    % of contours specified by Gamma.
    % Gamma is a 3D array of contours, parameterised by a set of point.
    % Each slice of Gamma is a separate contour, whose points x and y
    % coordinates are given by the first and second column.
    % The output cross_list is a cell array containing an array of all the
    % crossings for each point in x.

    Nc      = size(Gamma, 3); % number of contours
    Nx      = length(x);
    
    % Check if contours are closed, i.e. first and last points of each
    % contour is the same. If contours are not closed, this function will
    % not function properly for crossings around the first point.
    % Close contours if they are not already.
    if ~isequal(Gamma(1,:,:), Gamma(end,:,:))
        Gamma(end+1,:,:) = Gamma(1,:,:);
    end
    
    cross_list = cell(Nx, 1);
    for i = 1:Nx % for each sample point
        y_cross = [];

        for l = 1:Nc % for each contour
            X           = Gamma(:,1,l); % list of x-corrdinates in contour
            Y           = Gamma(:,2,l); % list of y-corrdinates in contour

            % Find all indices k for which x(i) is between X(k) and X(k+1)
            inside_inc  = (X <= x(i) & circshift(X,-1) > x(i)); % increasing X: X(k) <= x(i) < X(k+1)
            inside_dec  = (X >= x(i) & circshift(X,-1) < x(i)); % decreasing X: X(k) > x(i) > X(k+1)

            is_between  = (inside_inc | inside_dec);
            is_equal    = (X == x(i)); 

            % Remove points where the line goes straight through (these
            % give problems with the interpolation). Put them back later.
            to_interp   = is_between & ~is_equal;

            % Get the points in X and R "before" x(i)
            X0          = X(to_interp);
            Y0          = Y(to_interp);

            % Get the points in X and R "after" x(i)
            X1          = X(circshift(to_interp,1));
            Y1          = Y(circshift(to_interp,1));

            % Use linear interpolation to find the rapidity at all the
            % intersections between the contour and a vertical line going 
            % through x(i)
            y_cross     = [y_cross; ( Y0.*(X1-x(i)) + Y1.*(x(i)-X0) )./( X1 - X0 )];

            % Add back the points where the line goes straight through
            % Edge case! If the equal point is the first point of the
            % contour, it will be counted twice, as the contour is closed!
            is_equal(end) = false; % avoids the double counting
            y_cross    = [y_cross; Y(is_equal)];

            % If the line goes straight through a point and is tangent to
            % the countour, the point must be counted twice.
            % To find tangent point, note that its neighbours in x will be
            % on the same side of it.
            sign_left   = sign(X(1:end-1) - circshift(X(1:end-1),-1)); % must exclude last point as it is a duplicate of first, given contour is closed
            sign_right  = sign(X(1:end-1) - circshift(X(1:end-1),1));
            is_tangent  = sign_left == sign_right;
            Y           = Y(1:end-1); % remove last point to fit length of is_tangent            
            y_cross     = [y_cross; Y(is_tangent & is_equal(1:end-1))];
        end

        % sort the crossings from lowest to highest
        y_cross = sort(y_cross);

        % make sure that there is an equal number of crossings
        assert(mod(length(y_cross),2) == 0, 'Odd number of crossings!')

        % reshape to have sets of crossing, i.e. size(y_cross) = [2,N]
        y_cross = reshape(y_cross, 2, []);

        % Store the point x_ij and all its corresponding crossings
        cross_list{i} = y_cross;
    end
end